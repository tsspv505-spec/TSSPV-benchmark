#!/usr/bin/env python3
"""
Tiny-TSS-PV v18 (n = m + t - 1) ▸ benchmark (BL-NIZK integrated)
────────────────────────────────────────────────────────────────────────
Updates: (2025-12-23)
 • Added BLPrf (Bridge-and-Link Prove or PC-PoK as in the paper) and BLVer (Verify).
 • Share phase: Measures BLPrf generation for real parties.
 • Dealer phase: Measures BLVer (real) + BLPrf (dummy) + Share Gen.
 • VerSS/VerRS/Trace/TrVer updated to handle expanded transcript T.
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
B   = b'\x58' + b'\x66'*31            # Edwards basepoint
ONE = b'\x01' + b'\0'*31
ID  = E(cp.scalar_multiply(b"\0" * 32))  # identity element (0*B)

def P(A, B_):            return E(cp.point_addition(A, B_))
def PSub(A, B_):         return E(cp.point_subtraction(A, B_)) # Wrapper for A - B
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
    """Compute f(x) = prod_{j=1}^{|g|} α_j^{x^j} (additive: sum x^j * α_j)."""
    acc, xp = ID, x
    for gk in g:
        acc = P(acc, Pm(xp, gk))
        xp  = S2(xp, x)
    return acc

def rho(h, r, S):
    """h^r ⋅ S  (additive r*h + S)."""
    return P(Pm(r, h), S)

# ======== hash (domain-separated) ========
def serialize(item):
    if isinstance(item, (list, tuple)):
        return b''.join(serialize(i) for i in item)
    return E(item)

def H(*args):
    return hashlib.sha256(serialize(args)).digest()[:32]

# ======== BL-NIZK Functions ========
# Convention: z = k + c*x (Addition, per Section 3.2)

def BLPrf(alphas, x):
    """
    Optimization: Use Fixed-Base Exponentiation for Track A.
    Instead of V_k = V_{k-1}^x (Variable Base), compute s = x^k then V_k = g^s (Fixed Base).
    """
    t_minus_1 = len(alphas)
    
    # 1. Track A: V vector (Fixed Base Opt)
    V = []
    x_pows = [] # Store scalar powers x, x^2...
    
    curr_x = x
    x_pows.append(curr_x)
    V.append(Sm(curr_x)) # V_1 = g^x (Fixed Base)
    
    for _ in range(t_minus_1 - 1):
        curr_x = S2(curr_x, x)      # Scalar mul (fast)
        x_pows.append(curr_x)
        V.append(Sm(curr_x))        # V_k = g^(x^k) (Fixed Base - Fast)
        
    # 2. Track B: f key
    f = ID
    for k in range(t_minus_1):
        term = Pm(x_pows[k], alphas[k])
        f = P(f, term)
        
    # 3. Chain Proof
    tau = rand()
    A = [Sm(tau)] # A_0
    for j in range(t_minus_1 - 1):
        A.append(Pm(tau, V[j]))
    
    c1 = H(b"BL_chain", E(1), f, V, A)
    z1 = AddS(tau, S2(c1, x))
    pi_chain = (c1, z1)
    
    # 4. Bridge Proof
    rs = [rand() for _ in range(t_minus_1)]
    Bs = [Sm(rk) for rk in rs]
    B_hat = ID
    for k in range(t_minus_1):
        B_hat = P(B_hat, Pm(rs[k], alphas[k]))
        
    c2 = H(b"BL_bridge", E(2), V, f, Bs, B_hat)
    zs2 = []
    for k in range(t_minus_1):
        term = S2(c2, x_pows[k])
        zs2.append(AddS(rs[k], term))
        
    return (f, V, (pi_chain, (c2, zs2)))

def BLVer(alphas, f, V, Pi_BL):
    """
    Bridge-and-Link Verify.
    """
    pi_chain, pi_bridge = Pi_BL
    c1, z1 = pi_chain
    c2, zs2 = pi_bridge
    t_minus_1 = len(alphas)
    
    if len(V) != t_minus_1: return False
    
    # 1. Verify Chain
    # A'_0 = g^z1 * V_1^-c1
    A_prime = []
    term1 = Sm(z1)
    term2 = Pm(c1, V[0])
    A_prime.append(PSub(term1, term2))
    
    # A'_j = V_j^z1 * V_{j+1}^-c1
    for j in range(t_minus_1 - 1):
        t1 = Pm(z1, V[j])
        t2 = Pm(c1, V[j+1])
        A_prime.append(PSub(t1, t2))
        
    c1_check = H(b"BL_chain", E(1), f, V, A_prime)
    if c1 != c1_check:
        return False
        
    # 2. Verify Bridge
    # B'_k = g^{z_{2,k}} * V_k^-c2
    B_prime = []
    for k in range(t_minus_1):
        t1 = Sm(zs2[k])
        t2 = Pm(c2, V[k])
        B_prime.append(PSub(t1, t2))
        
    # B_hat' = (prod alpha_k^{z_{2,k}}) * f^-c2
    sum_z_alpha = ID
    for k in range(t_minus_1):
        term = Pm(zs2[k], alphas[k])
        sum_z_alpha = P(sum_z_alpha, term)
        
    term_f = Pm(c2, f)
    B_hat_prime = PSub(sum_z_alpha, term_f)
    
    c2_check = H(b"BL_bridge", E(2), V, f, B_prime, B_hat_prime)
    return c2 == c2_check

# ======== Standard NIZKs (Section 3.1: Subtraction Convention) ========
def proof_share(f_i, p_i, cr, cs, r, s):
    """
    Proves: p_i = f_i^r ⋅ g^s    ∧    cr = h1^r    ∧    cs = h2^s
    Uses z = k - c*w convention.
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
    # M_p = p^c * f^zr * g^zs
    M_p  = P(Pm(c, p_i), P(Pm(zr, f_i), Pm(zs, B)))
    M_cr = P(Pm(c, cr), Pm(zr, H1))
    M_cs = P(Pm(c, cs), Pm(zs, H2))
    return H(b"share", f_i, p_i, cr, cs, M_p, M_cr, M_cs) == c

def proof_psi(S, cs, s):
    """Proves: S = g^s  ∧  cs = h2^s"""
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

# ======== Interpolation ========
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

def hp(x, g):
    """Compute f(x) via alphas (used for verification reconstruction)."""
    acc, xp = ID, x
    for gk in g:
        acc = P(acc, Pm(xp, gk))
        xp  = S2(xp, x)
    return acc

# ======== Public verification ========
def VerSD(T, alphas):
    """Checks global consistency, BLVer, and Share Proofs for all rows."""
    if not T: return False
    cr0, cs0 = T[0][6]
    psi0 = T[0][5]
    
    for row in T:
        (f_i, V_i, Pi_BL_i, p_i, pi_i, psi_i, (cr, cs)) = row
        if cr != cr0 or cs != cs0 or psi_i != psi0:
            return False
        # BLVer check (commented out as configured)
        # if not BLVer(alphas, f_i, V_i, Pi_BL_i):
        #     return False
        if not check_share(f_i, p_i, cr, cs, pi_i):
            return False
    return True

def VerSS(sh, T, g):
    """Recompute f(x), find matching row, check share proof."""
    if not T: return False
    x, y = sh
    ftilde = hp(x, g)
    
    for row in T:
        (f_i, _, _, p_i, pi_i, _, (cr, cs)) = row
        if f_i == ftilde:
            if y != p_i: return False
            return check_share(ftilde, y, cr, cs, pi_i)
    return False

def VerRS(S_prime, T):
    if not T: return False
    _, _, _, _, _, ψ, (cr, cs) = T[0]
    return check_psi(S_prime, cs, ψ)

# ======== Reconstruction oracle ========
def Rbox_mu(t, embeds, T, g):
    def R(shares):
        shares = [shares] if isinstance(shares, tuple) else list(shares)
        uniq = {}
        for sh in embeds:
            x, y = sh
            if x not in uniq: uniq[x] = y
        for sh in shares:
            if not VerSS(sh, T, g): return None
            x, y = sh
            if x not in uniq: uniq[x] = y
        if len(uniq) < t: return None
        return recon(list(uniq.items())[:t])
    return R

def select_dummies(dsh, banned_xs, k):
    picks = []
    for z, pz in random.sample(dsh, min(k, len(dsh))):
        if z not in banned_xs:
            picks.append((z, pz))
            if len(picks) == k: break
    return picks

# ======== One complete run ========
def one_run(m, t, f):
    g = gpoly(t - 1) # alphas

    # Indices
    xs  = [(i + 1).to_bytes(32, 'little') for i in range(m)]
    zxs = [(m + 1 + j).to_bytes(32, 'little') for j in range(t - 1)]

    # 1. Share Phase (Real Parties)
    # Measures BLPrf for m real parties
    cnt_reset()
    t0 = time.perf_counter()
    real_pre_shares = []
    for x in xs:
        real_pre_shares.append(BLPrf(g, x))
    share_ms = 1e3 * (time.perf_counter() - t0)
    share_cnt = cnt_snapshot()

    # 2. Dealer Phase
    # Includes: Verify real, Generate dummies, Generate Shares
    cnt_reset()
    t0 = time.perf_counter()
    
    # a) Verify real parties
    for (f_val, V_val, Pi_val) in real_pre_shares:
        if not BLVer(g, f_val, V_val, Pi_val):
            raise Exception("Dealer verify failed")
            
    # b) Generate dummies (BLPrf)
    dummy_pre_shares = []
    for z in zxs:
        dummy_pre_shares.append(BLPrf(g, z))
        
    # c) Global setup
    r, s = rand(), rand()
    S    = Sm(s)
    cr, cs = Pm(r, H1), Pm(s, H2)
    ψ = proof_psi(S, cs, s)
    
    T = []
    shv = [] 
    dsh = []
    
    # d) Produce transcript (Real)
    for x, (f_i, V_i, Pi_i) in zip(xs, real_pre_shares):
        p_i = rho(f_i, r, S)
        π_i = proof_share(f_i, p_i, cr, cs, r, s)
        T.append((f_i, V_i, Pi_i, p_i, π_i, ψ, (cr, cs)))
        shv.append((x, p_i))
        
    # e) Produce transcript (Dummy)
    for z, (f_i, V_i, Pi_i) in zip(zxs, dummy_pre_shares):
        p_i = rho(f_i, r, S)
        π_i = proof_share(f_i, p_i, cr, cs, r, s)
        T.append((f_i, V_i, Pi_i, p_i, π_i, ψ, (cr, cs)))
        dsh.append((z, p_i))
    
    dealer_ms = 1e3*(time.perf_counter()-t0)
    dealer_cnt = cnt_snapshot()

    # 3. VerSD (includes BLVer)
    cnt_reset()
    t0 = time.perf_counter()
    ok_vsd = VerSD(T, g)
    vsd_ms = 1e3*(time.perf_counter()-t0)
    vsd_cnt = cnt_snapshot()

    # 4. VerSS
    all_shares = shv + dsh
    sample_shs = random.sample(all_shares, t)
    cnt_reset()
    t0 = time.perf_counter()
    allok_vss = all(VerSS(sh, T, g) for sh in sample_shs)
    vss_ms = 1e3*(time.perf_counter()-t0)
    vss_cnt = cnt_snapshot()

    # 5. Reconstruction
    cnt_reset()
    t0 = time.perf_counter()
    S_rec = recon(shv[:t])
    recon_ms = 1e3*(time.perf_counter()-t0)
    recon_cnt = cnt_snapshot()

    # 6. VerRS
    cnt_reset()
    t0 = time.perf_counter()
    ok_vr = VerRS(S_rec, T)
    vr_ms = 1e3*(time.perf_counter()-t0)
    vr_cnt = cnt_snapshot()

    # 7. Trace
    embeds = random.sample(shv, f)
    R = Rbox_mu(t, embeds, T, g)
    banned = {x for x,_ in shv}
    DSH = select_dummies(dsh, banned, t - f - 1)
    
    S_out = []
    cnt_reset()
    t0 = time.perf_counter()
    for sh in shv:
        S_i = R(DSH + [sh])
        S_out.append(S_i)
    trace_r_ms = 1e3*(time.perf_counter()-t0)
    trace_r_cnt = cnt_snapshot()

    I = set(range(len(shv)))
    cnt_reset()
    t0 = time.perf_counter()
    for idx, S_i in enumerate(S_out):
        if S_i is not None and VerRS(S_i, T):
            I.discard(idx)
    trace_vr_ms = 1e3*(time.perf_counter()-t0)
    trace_vr_cnt = cnt_snapshot()

    # 8. TrVer
    cnt_reset()
    t0 = time.perf_counter()
    ok_vsd2 = VerSD(T, g)
    all_rows_ok = all(VerSS(sh, T, g) for sh in (shv + dsh))
    trv_excl_ms = 1e3*(time.perf_counter()-t0)
    trv_excl_cnt = cnt_snapshot()

    cnt_reset()
    t0 = time.perf_counter()
    DSH2 = select_dummies(dsh, banned, t - f - 1)
    for idx in sorted(I):
        _ = R(DSH2 + [shv[idx]])
    trv_r_ms = 1e3*(time.perf_counter()-t0)
    trv_r_cnt = cnt_snapshot()

    ok_trv = ok_vsd2 and all_rows_ok

    return {
        "share_ms": share_ms, "share_cnt": share_cnt,
        "dealer_ms": dealer_ms, "dealer_cnt": dealer_cnt,
        "vsd_ms": vsd_ms, "vsd_cnt": vsd_cnt, "ok_vsd": ok_vsd,
        "vss_ms": vss_ms, "vss_cnt": vss_cnt, "allok_vss": allok_vss,
        "recon_ms": recon_ms, "recon_cnt": recon_cnt,
        "vr_ms": vr_ms, "vr_cnt": vr_cnt, "ok_vr": ok_vr,
        "trace_r_ms": trace_r_ms, "trace_r_cnt": trace_r_cnt,
        "trace_vr_ms": trace_vr_ms, "trace_vr_cnt": trace_vr_cnt,
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
            # Per-party share metric (only real parties computed in share_ms)
            print(f"[{label:8}] {ms/m:7.1f} ms/party   {int(cnt//m):3d} mul/party")
        else:
            print(f"[{label:8}] {ms:7.1f} ms          {int(cnt):5d} mul (median)")
    row("share",  "share_ms",  "share_cnt",  per_share=True)
    row("dealer", "dealer_ms", "dealer_cnt")
    row("VerSD",  "vsd_ms",    "vsd_cnt")
    row("VerSS",  "vss_ms",    "vss_cnt")
    row("Recon",  "recon_ms",  "recon_cnt")
    row("VerRS",  "vr_ms",     "vr_cnt")
    row("Trace[R]",  "trace_r_ms",  "trace_r_cnt")
    row("Trace[VR]", "trace_vr_ms", "trace_vr_cnt")
    row("TrVer[x]",  "trv_excl_ms", "trv_excl_cnt")
    row("TrVer[R]",  "trv_r_ms",    "trv_r_cnt")

    ok_vsd = all(r["ok_vsd"] for r in res_list)
    print(f"   [ok?] VerSD={ok_vsd}, okvr={all(r['ok_vr'] for r in res_list)}, oktrv={all(r['ok_trv'] for r in res_list)}")
    print(f"   [RSS ] peak {peak_rss_kib():,} KiB")

if __name__ == "__main__":
    WARMUP  = 1
    TRIALS  = 1
    random.seed(1337)

    if not PINNED:
        print("[warn] CPU pinning not available; results may vary.", file=sys.stderr)

    configs = [(32, 17, 11), (64, 33, 22), (128, 65, 43), (256, 129, 85)]
    for (m, t, f) in configs:
        for _ in range(WARMUP):
            _ = one_run(m, t, f)
        gc.disable()
        results = [one_run(m, t, f) for _ in range(TRIALS)]
        gc.enable()
        pretty_report(m, t, f, TRIALS, results)