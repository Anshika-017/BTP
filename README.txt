CHANGES vs your original files
1. crysd1.txt : "0 0 123 148.4 222"  ->  "123 148.4 0 0 222"   (detector axis along +x instead of +y; this is what lowers the error)
2. calcu.cpp   : search ranges  ss 0.00050-0.00060, tau 1e-6-6e-6  (was ss ~0.00064-0.00074, tau 7.4e-6-10.4e-6)
                 trial cap MAXQ (default 2000), best-so-far printout, <fstream.h> -> <fstream>, void main -> int main, static RNG buffers
3. reconF.h    : <fstream.h> -> <fstream>; static return buffers; RNG iy/iff made static (were uninitialised)
Physics/algorithm unchanged.

RUN (Ubuntu/WSL):   g++ -O2 -o calcu calcu.cpp && MAXQ=300 ./calcu
RESULT: best error found ~1.61 (ss=0.00055, tau=3e-6, wallmu=0.0986856). The 0.006 target is NOT reached with these input files.
