#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


int read_graph(const char* filename)
{
    int node1, node2, num_nodes = 0;
    float edge;
    FILE *pf;
    pf = fopen (filename, "r");
    if (pf == NULL)
        return 1;
    while (fscanf(pf, "%d %d %f", &node1, &node2, &edge) == 3)
    {
        num_nodes = MAX(node1, num_nodes);
        num_nodes = MAX(node2, num_nodes);
    }
    
    fclose (pf);
    printf("the number of nodes are %d \n", num_nodes);
    return 0;
}

int main(int argc, char *argv[])
{
    read_graph("./data/Graph_directed.txt");
    return 0;
}