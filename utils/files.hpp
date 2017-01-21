#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>





class File
{
private:
    
public:
    char name[256];
    char access_mode[8];
    FILE* stream;
    
    File(const char*, const char*);
    ~File();
};

File::File(const char* _name, const char* _access_mode)
{
    strcpy(name,_name);
    strcpy(access_mode,_access_mode);
    
    stream = fopen(name,access_mode);
    if (!stream)                        {  fprintf(stderr,"Error! Cannot open file "); exit(EXIT_FAILURE);  }
}


File::~File()
{
    if (!stream) fclose(stream);
}