#include <stdio.h>
#include <stdlib.h>

void cshift(double *data, int nx, int ny, int axis, int direction) {
    int x, y;
    double b;

    if (axis == 0 && direction == 1) {
        for (y = 0; y != ny; y++) {
            b = *(data + y);
            for (x = 0; x != nx - 1; x++) {
                *(data + y + ny * x) = *(data + y + ny * (x + 1));
            }
            *(data + y + ny * x) = b;
        }
    }

    if (axis == 0 && direction == -1) {
        for (y = 0; y != ny; y++) {
            b = *(data + y + ny * (nx - 1));
            for (x = nx - 1; x >= 0; x--) {
                *(data + y + ny * (x + 1)) = *(data + y + ny * x);
            }
            *(data + y) = b;
        }
    }

    if (axis == 1 && direction == 1) {
        for (x = 0; x < nx; x++) {
            b = *(data + x * ny);
            for (y = 0; y != ny; y++) {
                *(data + y + ny * x) = *(data + y + 1 + ny * x);
            }
            *(data + y - 1 + ny * x) = b;
        }
    }

    if (axis == 1 && direction == -1) {
        for (x = 0; x < nx; x++) {
            b = *(data + ny - 1 + ny * x);
            for (y = ny - 2; y >= 0; y--) {
                *(data + y + 1 + ny * x) = *(data + y + ny * x);
            }
            *(data + x * ny) = b;
        }
    }
}

void print_array(double *data, int nx, int ny) {
    int x, y;
    for (y = 0; y != ny; y++) {
        for (x = 0; x != nx; x++) {
            printf("%3i:%9.6f ",(int)(x * ny + y),data[x * ny + y]);
        }
        printf("\n");
    }
}

#include <limits.h>            /* for CHAR_BIT */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BITMASK(b) (1 << ((b) % CHAR_BIT))
#define BITSLOT(b) ((b) / CHAR_BIT)
#define BITSET(a, b) ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCLEAR(a, b) ((a)[BITSLOT(b)] &= ~BITMASK(b))
#define BITTEST(a, b) ((a)[BITSLOT(b)] & BITMASK(b))
#define BITNSLOTS(nb) ((nb + CHAR_BIT - 1) / CHAR_BIT)

void Print2DArray(double *A, int nr, int nc) {
    int r, c;
    for (r = 0; r < nr; r++) {
        for (c = 0; c < nc; c++)
            printf("%7.1f", *(A + r * nc + c));

        printf("\n");
    }
    printf("\n\n");
}

// Non-square matrix transpose of matrix of size r x c and base address data

void MatrixInplaceTranspose(double *data, int r, int c) {
    int size = r * c - 1;
    double temp1, temp2; // holds element to be replaced, eventually becomes next element to move
    int next; // location of 'temp1' to be moved
    int cycleBegin; // holds start of cycle
    int i; // iterator

    char bitarray[BITNSLOTS(size)];
    memset(bitarray, 0, BITNSLOTS(size));

    BITSET(bitarray, 0);
    BITSET(bitarray, size);

    i = 1; // Note that data[0] and data[size-1] won't move
    while (i < size) {
        cycleBegin = i;
        temp1 = *(data + i);
        do {
            next = (i * r) % size;
            //swap
            temp2 = *(data + next);
            *(data + next) = temp1;
            temp1 = temp2;
            BITSET(bitarray, i);
            i = next;
        } while (i != cycleBegin);

        // Get Next Move
        for (i = 1; i < size && BITTEST(bitarray, i); i++);
    }
}

// Driver program to test above function

int test_matrix(void) {
    int r = 5, c = 6;
    int size = r*c;
    double *A;
    int i;
    A = (double *) malloc(size * sizeof (double));

    for (i = 0; i < size; i++)
        A[i] = i + 11.1;

    Print2DArray(A, r, c);
    MatrixInplaceTranspose(A, r, c);
    Print2DArray(A, c, r);

    return 0;
}
