#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

// compile: gcc par_mat_mult.c -o par_mat_mult -lm -lpthread

// for time measuring in nano seconds
#define BILLION 1000000000L

// the max num of digits of a matrix element
#define MAXLEN 5

// the total number of file scanners: multiply two matrices only
#define SCANNERS 2

// the total number of multiplication workers: utilize parallism
#define WORKERS 8

typedef struct matrix_element {
    int x;
    int y;
    float val;
} item;

// arguments per scanner
typedef struct scan_dims_arguments {
    char* filename;
    int* x;
    int* y;
} scan_args;

// arguments per reader
typedef struct read_mat_arguments {
    char* filename;
    int transpose;
    item** mat;
} read_args;

// arguments per multiplier
typedef struct submult_mat_arguments {
    int i;
    item* lrow;
    int j;
    item* rcol;
    int len;
    item** new_mat;
} submult_args;

// arguments per mult worker
typedef struct mult_work_arguments {
    int qlen;
    submult_args* task_queue;
} worker_args;

// a first scan to determine the dimensions of a matrix
// orig prototype:
// void scan_dims(const char* filename, int* x, int* y);
void* scan_dims(void* data) {
    scan_args* mydata; 
    mydata = (scan_args*) data;

    const char* filename = mydata->filename;
    int* x = mydata->x;
    int* y = mydata->y;

    FILE* fptr = fopen(filename, "r");
    if(fptr == NULL) {
        perror("File open failed");
        exit(EXIT_FAILURE);
    }

    int rows = 0;
    // since the last item in a row is not followed by a comma
    int columns = 1;

    int fch;
    while ((fch = fgetc(fptr)) != EOF) {
        //putchar(fch);
        if (rows == 0) {
            if (fch == ',') {
                columns += 1;
            }
        }
        if (fch == '\n') {
            rows += 1;
        }
    }

    *x = rows;
    *y = columns;

    if (fclose(fptr) == EOF)
    {
        perror("File closing failed");
        exit(EXIT_FAILURE);
    }
    //printf("%s is a %dx%d matrix.\n", filename, rows, columns);
    pthread_exit((void*) data);
}


// convert string to int value
int conv_int(char digits[], int len) {
    int res = 0;
    int i;
    for (i=0; i<len; i++) {
        res += (digits[i]-'0')*pow( 10, (len-1-i) );
    }
    return res;
}

// convert string to float value
float conv_float(char digits[], int len) {
    float res = 0;
    int i;
    for (i=0; i<len; i++) {
        res += (digits[i]-'0')/pow( 10, i+1 );
    }
    return res;
}

// convert an integer to string 
void conv_char(float val, char* str, int maxlen, int* reallen) {
    int i, j;
    int start = 0;
    for (i=maxlen-1, j=0; i>=0; i--) {
        int base = pow(10, i);
        //printf("Testing base=%d\n", base);
        if ((!start) && (base > val)) {
            continue;
        }

        start = 1;
        int digit = (int)(val/base);
        val -= digit*base;
        str[j] = digit + '0'; 
        //printf("str[%d] %c\n", j, str[j]);
        j += 1;
    }
    *reallen = j;
    str[j]='\0';
}

// construct items in memory with known dims
// orig prototype:
// void read_mat(const char* filename, item* mat[], int transpose);
void* read_mat(void* data) {
    read_args* mydata;
    mydata = (read_args*) data;

    const char* filename = mydata->filename;
    item** mat = mydata->mat;
    int transpose = mydata->transpose; 

    FILE* fptr = fopen(filename, "r");
    if(fptr == NULL) {
        perror("File open failed");
        exit(EXIT_FAILURE);
    }

    int r = 0;
    int c = 0;

    int fch;
    int len = 0;
    char digits[MAXLEN];
    while ((fch = fgetc(fptr)) != EOF) {
        //putchar(fch);
        if ((fch == ',') || (fch == '\n')) {
            int value = conv_int(digits, len);
            item new_item = {
                .x = c,
                .y = r,
                .val = (float) value,
            }; 
            //printf("[%d,%d] %d", c, r, value);
            if (transpose) {
                // put the new item at the rth cell at cth col
                mat[c][r] = new_item;
            } else {
                // put the new item at the cth cell at rth row
                mat[r][c] = new_item;
            }

            len = 0;
            if (fch == '\n') {
                c = 0;
                r += 1;
                //printf("\n");
            } else {
                c += 1;
            }
        } else {
            digits[len] = fch;
            len += 1;
        }
    }

    if (fclose(fptr) == EOF)
    {
        perror("File closing failed");
        exit(EXIT_FAILURE);
    }
    pthread_exit((void*) data);
}

/* serial version
// task for each worker
float sub_mult(item* lrow, item* rcol, int len) {
    float new_val = 0;
    int i;
    for (i=0; i<len; i++) {
        new_val += (lrow[i].val)*(rcol[i].val);
        //printf("%d %d %.2f, ", lrow[i].x, lrow[i].y, lrow[i].val);
        //printf("%d %d %.2f\n", rcol[i].x, rcol[i].y, rcol[i].val);
    }
    //printf("%.2f\n", new_val);
    return new_val;
}

// left and new matrix is stored by rows, right matrix is stored by columns
void mats_mult(item* lmat[], item* rmat[], int x, int y, int l, item* new_mat[]) {
    int i,j;
    for (i=0; i<x; i++) {
        for (j=0; j<l; j++) {
            float new_val = sub_mult(lmat[i], rmat[j], y);
            item new_item = {
                .x = i,
                .y = j,
                .val = (int) new_val,
            };
            new_mat[i][j] = new_item;
            //printf("[%d,%d] %.1f\n", i, j, new_mat[i][j].val);
        }
    }
}
*/

// multiply one row with one column, put result in place
void sub_mult(int i, item* lrow, int j, item* rcol, int len, item* new_mat[]) {
    float new_val = 0;
    int l;
    for (l=0; l<len; l++) {
        new_val += (lrow[l].val)*(rcol[l].val);
    }
    
    item new_item = {
        .x = i,
        .y = j,
        .val = (int) new_val,
    };

    new_mat[i][j] = new_item;
}

// task for each worker
void* mult_work(void* data) {
    worker_args* mydata;
    mydata = (worker_args*) data;
    
    int qlen = mydata->qlen;
    submult_args* task_queue = mydata->task_queue;

    int t;
    for (t=0; t<qlen; t++) {
        int i = task_queue[t].i;
        item* lrow = task_queue[t].lrow;
        int j = task_queue[t].j;
        item* rcol = task_queue[t].rcol;
        int len = task_queue[t].len;
        item** new_mat = task_queue[t].new_mat;

        sub_mult(i, lrow, j, rcol, len, new_mat);
    }
    pthread_exit((void*) data);
}

// left and new matrix is stored by rows, right matrix is stored by columns
void mats_mult(item* lmat[], item* rmat[], int x, int y, int l, item* new_mat[]) {
    int i, j, id;

    // keep track of how long is each worker queue
    int work_dist[WORKERS] = {0};

    for (i=0; i<x; i++) {
        for (j=0; j<l; j++) {
            // worker id
            id = (10*i + j) % WORKERS;
            work_dist[id] += 1;
        }
    }

    worker_args part_data[WORKERS];
    submult_args* task_data[WORKERS];

    for (i=0; i<WORKERS; i++) {
        // allocate args space for each worker
        task_data[i] = (submult_args*) malloc( sizeof(submult_args) * work_dist[i]);
        
        worker_args mydata = {
            .qlen = work_dist[i], 
            .task_queue = task_data[i],
        };

        part_data[i] = mydata;
    }

    // reset task count to 0
    for (i=0; i<WORKERS; i++) {
        work_dist[i] = 0;
    }
    // allocate data to each worker
    for (i=0; i<x; i++) {
        for (j=0; j<l; j++) {
            // worker id
            id = (10*i + j) % WORKERS;
            submult_args subdata = {
                .i = i,
                .lrow = lmat[i],
                .j = j,
                .rcol = rmat[j],
                .len = y,
                .new_mat = new_mat,
            };

            // get current queue len
            int tail = work_dist[id];
            part_data[id].task_queue[tail] = subdata;
            work_dist[id] += 1;
        }
    }

    int rc;
    void* status;

    pthread_t workers[WORKERS];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for (i=0; i<WORKERS; i++) {
        rc = pthread_create(&workers[i], &attr, mult_work, (void*) &part_data[i]);
        if (rc) {
            printf("pthread_create failed with error code %d.\n", rc);
            exit (EXIT_FAILURE);
        }
    }

    pthread_attr_destroy(&attr);
    for (i=0; i<WORKERS; i++) {
        rc = pthread_join(workers[i], &status);
        if (rc) {
            printf("pthread_join failed with error code %d.\n", rc);
            exit (EXIT_FAILURE);
        }
    }

    for (i=0; i<WORKERS; i++) {
        free (task_data[i]);
    }

}

// write output to a file, new matrix is stored by rows
void write_mat(item* mat[], int x, int y, const char* filename) {
    FILE* fptr = fopen(filename, "w+");
    if (fptr == NULL) {
        perror("File open failed");
        exit(EXIT_FAILURE);
    }

    int i,j;
    char valstr[MAXLEN*2];
    for (i=0; i<x; i++) {
        for (j=0; j<y; j++) {
            float value; 
            value = mat[i][j].val;
            //printf("[%d,%d] %.1f\n", i, j, value);

            int len;
            conv_char(value, valstr, MAXLEN*2, &len);
            //printf("[%d,%d] %s\n", i, j, valstr);

            int c;
            for (c=0; c<len; c++) {
                putc(valstr[c], fptr);
            }
            if (j == y-1) {
                putc('\n', fptr);
            } else {
                putc(',', fptr);
            }
        }
    }
    
    if (fclose(fptr) == EOF)
    {
        perror("File closing failed");
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char* argv[]) {
    if (argc == 1) {
        printf("Usage: mat_mult left.mat right.mat result.mat\n");
        return EXIT_FAILURE;
    }

    int x, y, k, l;
    struct timespec start, end;

    scan_args scdata[SCANNERS];
    scan_args lmat_data = {
        .x = &x,
        .y = &y,
        .filename = argv[1],
    };
    
    scan_args rmat_data = {
        .x = &k,
        .y = &l,
        .filename = argv[2],
    };

    scdata[0] = lmat_data;
    scdata[1] = rmat_data;

    /* module 1: scan matrices */
    pthread_t scanners[SCANNERS]; 
    pthread_attr_t attr;
    
    int i, rc;
    void* status;
    clock_gettime(CLOCK_MONOTONIC, &start);

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    // spread out data to scanners
    for (i=0; i<SCANNERS; i++) {
        rc = pthread_create(&scanners[i], &attr, scan_dims, (void*) &scdata[i]);
        if (rc) {
            printf ("pthread_create failed with error code %d.\n", rc);
            return (EXIT_FAILURE);
        }
    }
    // free attribute and wait for other threads
    pthread_attr_destroy(&attr);
    for (i=0; i<SCANNERS; i++) {
        rc = pthread_join(scanners[i], &status);
        if (rc) {
            printf ("pthread_join failed with error code %d.\n", rc);
            return (EXIT_FAILURE);
        }
    }
    
    
    clock_gettime(CLOCK_MONOTONIC, &end);
    long time_elapsed = BILLION*(end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec);
    
    printf("input\t%d numbers\n", (x*y + k*l));
    printf("output\t%d numbers\n", (x*l));
    printf("------------\n");
    printf("scan:\t%ld ns\n", time_elapsed);
    

    if (y != k) {
        printf("Warning: cannot multiply matrices from %s and %s.\n", argv[1], argv[2] );
        return EXIT_FAILURE;
    }

    /* module 2: read matrices into memory*/
    // the left matrix is stored by rows
    item* lmat[x];
    // the right matrix is stored by columns
    item* rmat[l];

    int r;
    for (r=0; r<x; r++) {
        lmat[r] = (item*) malloc(sizeof(item)*y);
    }
    int c;
    for (c=0; c<l; c++) {
        rmat[c] = (item*) malloc(sizeof(item)*y);
    }

    read_args rddata[SCANNERS];
    read_args lmat_rddata = {
        .filename = argv[1],
        .transpose = 0,
        .mat = lmat,
    }; 
    
    read_args rmat_rddata = {
        .filename = argv[2],
        .transpose = 1,
        .mat = rmat,
    }; 

    rddata[0] = lmat_rddata;
    rddata[1] = rmat_rddata;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    // spread data to readers
    for (i=0; i<SCANNERS; i++) {
        rc = pthread_create(&scanners[i], &attr, read_mat, (void*) &rddata[i]);
        if (rc) {
            printf ("pthread_create failed with error code %d.\n", rc);
            return (EXIT_FAILURE);
        }
    }
    pthread_attr_destroy(&attr);
    for (i=0; i<SCANNERS; i++) {
        rc = pthread_join(scanners[i], &status);
        if (rc) {
            printf ("pthread_join failed with error code %d.\n", rc);
            return (EXIT_FAILURE);
        }
    }
    
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_elapsed = BILLION*(end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec);
    printf("read:\t%ld ns\n", time_elapsed);
    
    

    /* module 3: execute the multiplication operation */
    // the new matrix is stored by rows
    item* new_mat[x];
    for (r=0; r<x; r++) {
        new_mat[r] = (item*) malloc(sizeof(item)*l);

    }

    clock_gettime(CLOCK_MONOTONIC, &start);

    mats_mult(lmat, rmat, x, y, l, new_mat);

    clock_gettime(CLOCK_MONOTONIC, &end);

    time_elapsed = BILLION*(end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec);
    printf("mult:\t%ld ns\n", time_elapsed);
    
    /* module 4: write out the result */
    clock_gettime(CLOCK_MONOTONIC, &start);

    write_mat(new_mat, x, l, argv[3]);

    clock_gettime(CLOCK_MONOTONIC, &end);

    time_elapsed = BILLION*(end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec);
    printf("write:\t%ld ns\n", time_elapsed);
    
    // clean up
    for (r=0; r<x; r++) {
        free (lmat[r]);
        free (new_mat[r]);
    }
    for (c=0; c<l; c++) {
        free (rmat[c]);
    }
    
    pthread_exit(NULL);
}
