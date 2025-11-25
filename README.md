# Traceable Secret Sharing with Public Verifiability (TSS-PV) Benchmark
1. Prerequisites
   - AVX2 enabled CPU
   - Python 3.10
   - Rustc 1.83.0
   - Cargo 1.83.0
   - Maturin 1.83
   - Curve25519-dalek
2. Run the build script 
 - IMPORTANT NOTE: We strongly recommend running TSS-PV with AVX2 enabled. The build script assumes that AVX2 is enabled to match the benchmark results in TSS-PV paper.
 - Without AVX2 enabled, remove `RUSTFLAGS='--cfg curve25519_dalek_backend="simd" -C target_feature=+avx2' ` from the `./build_script.sh` before running (NOTE: without AVX2, the performance will be largely different with the benchmark reported in the paper.)
    ```
    cd ./curve25519_python
    ./build_script.sh
    ``` 
3. Run the benchmark: `python3 ./benchmark/benchmark.py`
   
