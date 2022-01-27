#ifndef FASTIO_HPP
#define FASTIO_HPP

#include <cstdio>
#include <string>
#include "constants.hpp"
#include "type.hpp"
#include "filesystem.hpp"

template<class Tp>
class fastIO {
private:
    Tp *buffer, *S = nullptr, *T = nullptr;
    FILE * f;
	long bytes = 0;
    int sz;

public: 
    fastIO(FILE * f_) {
        f = f_;
        buffer = new Tp[PAGESIZE];
        sz = sizeof(Tp);
    }

    fastIO(const std::string &file, const std::string openStyle) {
        f = fopen(file.c_str(), openStyle.c_str());
        buffer = new Tp[PAGESIZE];
		if(openStyle[0] == 'r') {
			bytes = file_size(file);
		}
        sz = sizeof(Tp);
    }

    ~fastIO() {
        fclose(f);
        delete [] buffer;
    }

    char getChar()  
    {  
        if(S==T)  
        {   
			if(bytes < PAGESIZE) {
				T=(S=buffer) + fread(buffer,1, bytes, f);
			}
            else {
				T=(S=buffer) + fread(buffer,1,PAGESIZE, f);
				bytes -= PAGESIZE;
			}
            if(S==T) return EOF;  
        }
        return *S++;  
    }
    void getInt(int & re)
    {
        char c;   
        for(c=getChar();c<'0'||c>'9';c=getChar());  
        while(c>='0'&&c<='9')  
            re=(re<<1)+(re<<3)+(c-'0'),c=getChar();
    }
    void getUInt(unsigned & re)
    {
        char c;
        re = 0;
        for(c=getChar();c<'0'||c>'9';c=getChar());  
        while(c>='0'&&c<='9')  
            re=(re<<1)+(re<<3)+(c-'0'),c=getChar();
    }

    Tp getFromBin()
    {
		if(S==T)  
        {
            if(bytes < PAGESIZE * sz) {
				T=(S=buffer) + fread(buffer, sz, bytes/sz, f);
			}
            else {
				T=(S=buffer) + fread(buffer, sz, PAGESIZE, f);
				bytes -= PAGESIZE * sz;
			}
        }
        return *S++;
    }

    void writeArray(Tp * v, v_size n) {
        fwrite(v, sizeof(Tp), n, f);
    }

    void seek(v_size b) {
        assert(fseek(f, b, SEEK_SET) == 0);
    }

    void getFromBinRandom(Tp * p, v_size n, v_size offset) {
        assert(fseek(f, offset, SEEK_SET) == 0);
        assert(fread(p, 4, n, f) == n);
    }

    // void close() {
    //     fclose(f);
    // }
};

#endif