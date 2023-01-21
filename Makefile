CC=gcc

CFLAGS=-march=native -g3 -mavx2 -funroll-loops

OFLAGS=-Ofast -fopenmp 
# -fopt-info-all=nbody.gcc.optrpt



all: nbody3D nbody3D[V1] nbody3D[V2] nbody3D[V3] nbody3D[V4] nbody3D[V5] nbody3D[V6]

nbody3D: nbody.c
	gcc $(CFLAGS) $(OFLAGS) $< -o $@gcc -lm
	clang $(CFLAGS) $(OFLAGS) $< -o $@clang -lm

nbody3D[V1]: [V1]nbody3D_SoA.c
	gcc $(CFLAGS) $(OFLAGS) $< -o $@gcc -lm
	clang $(CFLAGS) $(OFLAGS) $< -o $@clang -lm

nbody3D[V2]: [V2]nbody3D_noPow.c
	gcc $(CFLAGS) $(OFLAGS) $< -o $@gcc -lm
	clang $(CFLAGS) $(OFLAGS) $< -o $@clang -lm

nbody3D[V3]: [V3]nbody3D_memory_align.c
	gcc $(CFLAGS) $(OFLAGS) $< -o $@gcc -lm
	clang $(CFLAGS) $(OFLAGS) $< -o $@clang -lm

nbody3D[V4]: [V4]nbody3D_unroll.c
	gcc $(CFLAGS) $(OFLAGS) $< -o $@gcc -lm
	clang $(CFLAGS) $(OFLAGS) $< -o $@clang -lm

nbody3D[V5]: [V5]nbody3D_auto_vecto.c
	gcc $(CFLAGS)  $(OFLAGS) $< -o $@gcc -lm
	clang $(CFLAGS)  $(OFLAGS) $< -o $@clang -lm


nbody3D[V6]: [V6]nbody3D_avx2.c
	gcc $(CFLAGS) $(OFLAGS) $< -o $@gcc -lm
	clang $(CFLAGS) $(OFLAGS) $< -o $@clang -lm

run_gcc: nbody3D[V1]
	taskset -c 3 ./nbody3D[V1]gcc
	taskset -c 3 ./nbody3D[V2]gcc
	taskset -c 3 ./nbody3D[V3]gcc
	taskset -c 3 ./nbody3D[V4]gcc
	taskset -c 3 ./nbody3D[V5]gcc
	taskset -c 3 ./nbody3D[V6]gcc

run_clang: nbody3D[V1]
	taskset -c 3 ./nbody3D[V1]clang
	taskset -c 3 ./nbody3D[V2]clang
	taskset -c 3 ./nbody3D[V3]clang
	taskset -c 3 ./nbody3D[V4]clang
	taskset -c 3 ./nbody3D[V5]clang
	taskset -c 3 ./nbody3D[V6]clang

clean:
	rm -Rf *~ nbody3D *.optrpt


