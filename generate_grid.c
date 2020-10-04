#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define FILE_NAME "initial_grid.txt"
#define DEAD '0'
#define ALIVE '1'

int main(int argc, char **argv){
    FILE* fp = fopen(FILE_NAME, "w");

    if(argc < 2){
        printf("Grid width is required but none was given!\n");
        return 1;
    }

    srand(time(0));
    int width = atoi(argv[1]);
    for(int i = 0; i < width; i++){
        for(int j = 0; j < width; j++){
            if((rand()% 100) < 50)
                fputc(DEAD, fp);
            else
                fputc(ALIVE, fp);
        }
        fputc('\n', fp);
    }
    fclose(fp);
    return 0;
}