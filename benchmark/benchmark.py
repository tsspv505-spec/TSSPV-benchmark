#!/usr/bin/env python3
"""
Tiny-TSS-PV v16 (n = m + t - 1) ▸ benchmark (analysis-aligned counters)
────────────────────────────────────────────────────────────────────────
This revision aligns per-phase #mul with the paper’s Efficiency Analysis:
 • Recon counts only t mults; VerRS counted separately (4 mults)
 • Trace is split into R-only (box) and VerRS (public checks)
 • TrVer is split into excl-box (VerSD+VerSS) and R-only (box replays)
Other features preserved: medians, CPU pinning (best-effort), peak RSS.
"""

import os, sys, platform, secrets, hashlib, time, random, gc, statistics, resource
import curve25519_python as cp

# ======== (optional) CPU pinning ========
def pin_to_first_cpu():
    try:
        if hasattr(os, "sched_setaffinity"):
            cpus = os.sched_getaffinity(0)
            first = min(cpus)
            os.sched_setaffinity(0, {first})
            return True
        return False
    except Exception:
        return False

PINNED = pin_to_first_cpu()

# ======== process memory (peak RSS) ========
def peak_rss_kib():
    """Peak RSS in KiB (ru_maxrss is KiB on Linux, bytes on macOS/BSD)."""
    ru = resource.getrusage(resource.RUSAGE_SELF)
    val = ru.ru_maxrss
    if platform.system() == "Darwin":
        return val // 1024
    return val

# ======== deterministic enc helpers ========
def E(v):  # normalize to bytes
    return v if isinstance(v, (bytes, bytearray)) else bytes(v)

# ======== global counters ========
CNT = 0
def cnt_reset():
    global CNT
    CNT = 0
def cnt_snapshot():
    return CNT

# ======== Curve wrappers ========
B   = b'\x58' + b'\x66'*31           # Edwards basepoint
ONE = b'\x01' + b'\0'*31
ID  = E(cp.scalar_multiply(b"\0" * 32))  # identity element (0*B)

def P(A, B_):            return E(cp.point_addition(A, B_))
def Pm(s, P_):
    global CNT
    CNT += 1
    return E(cp.point_multiply(s, P_))
def Sm(s):
    global CNT
    CNT += 1
    return E(cp.scalar_multiply(s))
def S2(a, b):            return E(cp.scalar_multiply_scalar(a, b))
def Inv(z):              return E(cp.scalar_inverse(z))
def AddS(a, b):          return E(cp.scalar_addition(a, b))
def Sub(a, b):           return E(cp.scalar_subtraction(a, b))
def rand():              return secrets.token_bytes(32)

# Generators for commitments (counted outside measured windows)
H1, H2 = Pm(b'\x02'+b'\0'*31, B), Pm(b'\x03'+b'\0'*31, B)

# ======== polynomial & sharing helpers ========
def gpoly(t):
    """Return α_1..α_{t-1} as group elements (points)."""
    return [Sm(rand()) for _ in range(t)]

def hp(x, g):
    """Compute f(x) = prod_{j=1}^{|g|} α_j^{x^j} in additive group form."""
    acc, xp = ID, x
    for gk in g:
        acc = P(acc, Pm(xp, gk))
        xp  = S2(xp, x)
    return acc

def rho(h, r, S):
    """h^r ⋅ S  (additive r*h + S)."""
    return P(Pm(r, h), S)

# ======== hash (domain-separated) ========
def H(*elts):
    return hashlib.sha256(b''.join(elts)).digest()[:32]

# ======== NIZKs ========
def proof_share(f_i, p_i, cr, cs, r, s):
    """
    Proves: p_i = f_i^r ⋅ g^s   ∧   cr = h1^r   ∧   cs = h2^s
    FS challenge binds (f_i, p_i, cr, cs) with a domain tag.
    """
    kr, ks = rand(), rand()
    A_p  = P(Pm(kr, f_i), Pm(ks, B))
    A_cr = Pm(kr, H1)
    A_cs = Pm(ks, H2)
    c = H(b"share", f_i, p_i, cr, cs, A_p, A_cr, A_cs)
    zr = Sub(kr, S2(c, r))
    zs = Sub(ks, S2(c, s))
    return (c, zr, zs)

def check_share(f_i, p_i, cr, cs, π):
    c, zr, zs = π
    M_p  = P(Pm(c, p_i), P(Pm(zr, f_i), Pm(zs, B)))
    M_cr = P(Pm(c, cr), Pm(zr, H1))
    M_cs = P(Pm(c, cs), Pm(zs, H2))
    return H(b"share", f_i, p_i, cr, cs, M_p, M_cr, M_cs) == c

def proof_psi(S, cs, s):
    """Proves: S = g^s  ∧  cs = h2^s (global, counted outside dealer window)."""
    k = rand()
    A_S  = Pm(k, B)
    A_cs = Pm(k, H2)
    c = H(b"psi", S, cs, A_S, A_cs)
    z = Sub(k, S2(c, s))
    return (c, z)

def check_psi(S_prime, cs, ψ):
    c, z = ψ
    M_S  = P(Pm(c, S_prime), Pm(z, B))
    M_cs = P(Pm(c, cs),      Pm(z, H2))
    return H(b"psi", S_prime, cs, M_S, M_cs) == c

# ======== Interpolation (Δ_j at 0) ========
def montgomery_batch_invert(values):
    k = len(values)
    partials = [ONE]
    for v in values:
        partials.append(S2(partials[-1], v))
    inv_total = Inv(partials[-1])
    invs = [None] * k
    for i in range(k-1, -1, -1):
        invs[i] = S2(partials[i], inv_total)
        inv_total = S2(inv_total, values[i])
    return invs

def recon(shs):
    """
    shs: list[(x_j, y_j)]
    Computes S = ∏_j y_j^{Δ_j} with Δ_j = ∏_{i≠j} x_i / (x_i - x_j)
    """
    acc = ID
    k = len(shs)
    denoms = [ONE] * k
    for j in range(k):
        for m in range(k):
            if m != j:
                denoms[j] = S2(denoms[j], Sub(shs[m][0], shs[j][0]))
    inv_denoms = montgomery_batch_invert(denoms)
    for j, (xj, yj) in enumerate(shs):
        l = inv_denoms[j]
        for m in range(k):
            if m != j:
                l = S2(l, shs[m][0])
        acc = P(acc, Pm(l, yj))
    return acc

# ======== Public verification ========
def VerSD(T):
    """All rows share same (cr,cs) and each π_i verifies."""
    if not T: return False
    cr0, cs0 = T[0][4]
    for (f_i, p_i, π_i, ψ_i, (cr, cs)) in T:
        if cr != cr0 or cs != cs0:
            return False
        if not check_share(f_i, p_i, cr, cs, π_i):
            return False
    return True

def VerSS(sh, T, g):
    """Recompute f(x), find matching row, check y==p_i and verify π_i."""
    if not T: return False
    x, y = sh
    ftilde = hp(x, g)
    for (f_i, p_i, π_i, ψ_i, (cr, cs)) in T:
        if f_i == ftilde:
            if y != p_i:
                return False
            return check_share(ftilde, y, cr, cs, π_i)
    return False

def VerRS(S_prime, T):
    """Check global ψ on (S', cs)."""
    if not T: return False
    _, _, _, ψ, (cr, cs) = T[0]
    return check_psi(S_prime, cs, ψ)

# ======== Reconstruction oracle (Rbox) ========
def Rbox_mu(t, embeds, T, g):
    def R(shares):
        shares = [shares] if isinstance(shares, tuple) else list(shares)
        # Keep embedded first; prevent overriding by later duplicates
        uniq = {}
        for sh in embeds:
            x, y = sh
            if x not in uniq:
                uniq[x] = y
        # Validate submitted shares, build uniq set
        for sh in shares:
            if not VerSS(sh, T, g):
                return None
            x, y = sh
            if x not in uniq:
                uniq[x] = y
        if len(uniq) < t:
            return None
        return recon(list(uniq.items())[:t])
    return R

# ======== Dummy selection ========
def select_dummies(dsh, banned_xs, k):
    picks = []
    for z, pz in random.sample(dsh, min(k, len(dsh))):
        if z not in banned_xs:
            picks.append((z, pz))
            if len(picks) == k:
                break
    return picks

# ======== One complete run (returns per-phase timings & counts) ========
def one_run(m, t, f):
    # α_1..α_{t-1}
    g = gpoly(t - 1)

    # Party and dummy indices (pairwise distinct)
    xs  = [(i + 1).to_bytes(32, 'little') for i in range(m)]
    zxs = [(m + 1 + j).to_bytes(32, 'little') for j in range(t - 1)]

    # Parties/dummies compute f(x)
    cnt_reset()
    t0 = time.perf_counter()
    fx = [hp(x, g)  for x in xs]
    fz = [hp(z, g)  for z in zxs]
    share_ms = 1e3 * (time.perf_counter() - t0)
    share_cnt = cnt_snapshot()

    # Dealer picks r,s; global S, C, ψ (global cost left uncounted by design)
    r, s = rand(), rand()
    S    = Sm(s)                 # global
    cr, cs = Pm(r, H1), Pm(s, H2)  # global
    ψ = proof_psi(S, cs, s)      # global

    # Dealer per-row work: p_i and π_i
    cnt_reset()
    t0 = time.perf_counter()
    T = []
    for x, f_i in zip(xs, fx):  # parties
        p_i = rho(f_i, r, S)
        π_i = proof_share(f_i, p_i, cr, cs, r, s)
        T.append((f_i, p_i, π_i, ψ, (cr, cs)))
    dsh = []
    for z, f_i in zip(zxs, fz):  # dummies
        p_i = rho(f_i, r, S)
        π_i = proof_share(f_i, p_i, cr, cs, r, s)
        T.append((f_i, p_i, π_i, ψ, (cr, cs)))
        dsh.append((z, p_i))
    shv = [(x, p) for (x, (f, p, _, _, _)) in zip(xs, T)]  # (x, p)
    dealer_ms = 1e3*(time.perf_counter()-t0)
    dealer_cnt = cnt_snapshot()

    # VerSD (7n mults)
    cnt_reset()
    t0 = time.perf_counter()
    ok_vsd = VerSD(T)
    vsd_ms = 1e3*(time.perf_counter()-t0)
    vsd_cnt = cnt_snapshot()

    # VerSS on t random shares from all rows (t*(t+6) mults)
    all_shares = shv + dsh
    sample_shs = random.sample(all_shares, t)
    cnt_reset()
    t0 = time.perf_counter()
    allok_vss = all(VerSS(sh, T, g) for sh in sample_shs)
    vss_ms = 1e3*(time.perf_counter()-t0)
    vss_cnt = cnt_snapshot()

    # Reconstruction ONLY (t mults)
    cnt_reset()
    t0 = time.perf_counter()
    S_rec = recon(shv[:t])
    recon_ms = 1e3*(time.perf_counter()-t0)
    recon_cnt = cnt_snapshot()

    # VerRS ONLY (4 mults)
    cnt_reset()
    t0 = time.perf_counter()
    ok_vr = VerRS(S_rec, T)
    vr_ms = 1e3*(time.perf_counter()-t0)
    vr_cnt = cnt_snapshot()

    # R box and tracing — split into R-only vs VerRS(public)
    embeds = random.sample(shv, f)
    R = Rbox_mu(t, embeds, T, g)
    banned = {x for x,_ in shv}
    DSH = select_dummies(dsh, banned, t - f - 1)
    assert len(DSH) == t - f - 1, "insufficient dummy shares"

    # Trace: R-only
    S_out = []
    cnt_reset()
    t0 = time.perf_counter()
    for sh in shv:
        S_i = R(DSH + [sh])
        S_out.append(S_i)
    trace_r_ms = 1e3*(time.perf_counter()-t0)
    trace_r_cnt = cnt_snapshot()

    # Trace: VerRS over non-leaker outputs; produce I and Π
    I = set(range(len(shv)))
    cnt_reset()
    t0 = time.perf_counter()
    for idx, S_i in enumerate(S_out):
        if S_i is not None and VerRS(S_i, T):
            I.discard(idx)
    trace_vr_ms = 1e3*(time.perf_counter()-t0)
    trace_vr_cnt = cnt_snapshot()
    Π = [shv[i] for i in sorted(I)]

    # TrVer: excl-box = VerSD + VerSS(all rows)
    cnt_reset()
    t0 = time.perf_counter()
    ok_vsd2 = VerSD(T)
    all_rows_ok = all(VerSS(sh, T, g) for sh in (shv + dsh))
    trv_excl_ms = 1e3*(time.perf_counter()-t0)
    trv_excl_cnt = cnt_snapshot()
    ok_trv_prefix = ok_vsd2 and all_rows_ok

    # TrVer: R-only replays on accused (no VerRS if accusations are correct)
    cnt_reset()
    t0 = time.perf_counter()
    # Fresh DSH for replay
    DSH2 = select_dummies(dsh, banned, t - f - 1)
    for idx in sorted(I):
        _ = R(DSH2 + [shv[idx]])  # should return None for colluders
    trv_r_ms = 1e3*(time.perf_counter()-t0)
    trv_r_cnt = cnt_snapshot()

    ok_trv = ok_trv_prefix  # in the success case, R returns None on replays

    return {
        "share_ms": share_ms, "share_cnt": share_cnt,
        "dealer_ms": dealer_ms, "dealer_cnt": dealer_cnt,
        "vsd_ms": vsd_ms, "vsd_cnt": vsd_cnt, "ok_vsd": ok_vsd,
        "vss_ms": vss_ms, "vss_cnt": vss_cnt, "allok_vss": allok_vss,
        "recon_ms": recon_ms, "recon_cnt": recon_cnt,
        "vr_ms": vr_ms, "vr_cnt": vr_cnt, "ok_vr": ok_vr,
        "trace_r_ms": trace_r_ms, "trace_r_cnt": trace_r_cnt,
        "trace_vr_ms": trace_vr_ms, "trace_vr_cnt": trace_vr_cnt,
        "I_len": len(I),
        "trv_excl_ms": trv_excl_ms, "trv_excl_cnt": trv_excl_cnt,
        "trv_r_ms": trv_r_ms, "trv_r_cnt": trv_r_cnt, "ok_trv": ok_trv,
    }

# ======== Aggregate utilities ========
def med(values): return statistics.median(values)

def pretty_report(m, t, f, trials, res_list):
    n_total = m + (t - 1)
    print(f"\n== n={n_total} (total) | m={m} (real) | t={t} | f={f} | trials={trials} | pinned={PINNED} ==")
    def row(label, key_ms, key_cnt, per_share=False):
        ms = med([r[key_ms] for r in res_list])
        cnt = med([r[key_cnt] for r in res_list])
        if per_share:
            print(f"[{label:8}] {ms/n_total:7.1f} ms/share   {int(cnt//n_total):3d} mul/share")
        else:
            print(f"[{label:8}] {ms:7.1f} ms         {int(cnt):5d} mul (median)")
    # Sharing & dealer
    row("share",  "share_ms",  "share_cnt",  per_share=True)
    row("dealer", "dealer_ms", "dealer_cnt")
    # Public verification
    row("VerSD",  "vsd_ms",    "vsd_cnt")
    row("VerSS",  "vss_ms",    "vss_cnt")
    # Reconstruction vs VerRS (separated)
    row("Recon",  "recon_ms",  "recon_cnt")
    row("VerRS", "vr_ms",     "vr_cnt")
    # Trace: split
    row("Trace[R]",  "trace_r_ms",  "trace_r_cnt")
    row("Trace[VR]", "trace_vr_ms", "trace_vr_cnt")
    # TrVer: split
    row("TrVer[x]",  "trv_excl_ms", "trv_excl_cnt")
    row("TrVer[R]",  "trv_r_ms",    "trv_r_cnt")
    # Status
    ok_vsd = all(r["ok_vsd"] for r in res_list)
    ok_vr  = all(r["ok_vr"]  for r in res_list)
    ok_trv = all(r["ok_trv"] for r in res_list)
    print(f"   [ok?] VerSD={ok_vsd}  VerRS={ok_vr}  TrVer={ok_trv}")
    print(f"   [RSS ] peak {peak_rss_kib():,} KiB")

# ======== Main ========
if __name__ == "__main__":
    WARMUP  = 1
    TRIALS  = 10
    random.seed(1337)

    if not PINNED:
        print("[warn] CPU pinning not available; results may vary. Consider `taskset -c 0 python3 ...`.", file=sys.stderr)

    configs = [(32, 17, 11), (64, 33, 22), (128, 65, 43), (256, 129, 85)]
    for (m, t, f) in configs:
        for _ in range(WARMUP):
            _ = one_run(m, t, f)
        gc.disable()
        results = [one_run(m, t, f) for _ in range(TRIALS)]
        gc.enable()
        pretty_report(m, t, f, TRIALS, results)
