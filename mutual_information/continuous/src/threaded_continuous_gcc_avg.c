#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <pthread.h>

/* Constants */
#define MAX_NUM_VALUES 267600
#define MAX_VARIABLES 95
#define MAX_K_VALUE 5

/* Struct definitions */
typedef struct data_shape {
    int num_variables;
    int num_entries_per_variable;
} DataShape;

typedef struct indexed_value {
    double value;
    int index;
} IndexedValue;

typedef struct point {
    double coords[2];
} Point;

typedef struct kdtree {
    double value;
    int middle_index;
    int left_index;
    int right_index;
} KDTree;

typedef struct combine_mi_thread_args {
    int *valid_values;
    double *psi_values;
    double *data;
    IndexedValue *data_value_arrays;
    double *data_matrix;
    int num_attributes;
    int num_variables;
    int num_entries_per_variable;
    int k;
    int p;
} CombineMiThreadArgs;

/* Function type declarations */
void find_knearest_neighbour_distance_helper(KDTree *kdtree_array, Point *point_array, Point point, double *min_distance_array, int k, int axis, int lower_index, int higher_index, int curr_tree_index);

/* --- Private Implementation Functions (should not be used directly) --- */

/*
 * Comparison function to be used for sorting an array of IndexedValues in
 * ascending order using quick sort
 */
int compare_value(const void *a, const void *b) {
    double a_value = ((IndexedValue *) a)->value;
    double b_value = ((IndexedValue *) b)->value;
    return (a_value > b_value) - (a_value < b_value);
}

/*
 * Optimized implementation of the Digamma function restricted to postive
 * values of x.
 */
double psi(double x) {
    double r = 0;
    while (x <= 5) {
        r -= 1/x;
        x += 1;
    }

    double f = 1 / (x * x);
    double t = f*(-1/12.0 + f*(1/120.0 + f*(-1/252.0 + f*(1/240.0 + f*(-1/132.0 + f*(691/32760.0 + f*(-1/12.0 + f*3617/8160.0)))))));
    return r + log(x) - 0.5/x + t;
}

/*
 * Computes the max norm between two 2D points
 */
double max_norm(Point point1, Point point2) {
    double x_dist = fabs(point1.coords[0] - point2.coords[0]);
    double y_dist = fabs(point1.coords[1] - point2.coords[1]);
    return (x_dist > y_dist) ? x_dist : y_dist;
}

double query_ball_point_fast(IndexedValue *value_array, int length, double *radius_array, double *psi_values) {
    double total = 0;
    for (int i = 0; i < length; i++) {
        int step;
        int lower_index;
        int middle_index;
        int higher_index;
        double radius = radius_array[value_array[i].index];

        // Find lower index
        step = 1;
        lower_index = i;
        int left_index = -1;
        while (lower_index >= 0 && value_array[i].value - value_array[lower_index].value < radius) {
            lower_index -= step;
            step = step << 1;
        }

        // Reverse previous step to get higher index
        step = step >> 1;
        higher_index = lower_index + step;

        // Perform iterative binary search to find left index (left and right indices are inclusive)
        while (left_index == -1) {
            middle_index = (lower_index + higher_index) >> 1;
            // Check if middle index is valid
            if (middle_index >= 0 && value_array[i].value - value_array[middle_index].value < radius) {
                // Check if the index to the left of the middle index is invalid
                if ((middle_index - 1) < 0 || value_array[i].value - value_array[middle_index - 1].value >= radius) {
                    // Found
                    left_index = middle_index;
                } else {
                    // Search left
                    higher_index = middle_index - 1;
                }
            } else {
                // Search right
                lower_index = middle_index + 1;
            }
        }

        // Find starting boundary indices for binary search
        step = 1;
        higher_index = i;
        int right_index = -1;
        while (higher_index < length && value_array[higher_index].value - value_array[i].value < radius) {
            higher_index += step;
            step = step << 1;
        }

        // Reverse previous step to get higher index
        step = step >> 1;
        lower_index = higher_index - step;

        // Perform iterative binary search to find left index (left and right indices are inclusive)
        while (right_index == -1) {
            middle_index = (lower_index + higher_index) >> 1;
            // Check if middle index is valid
            if (middle_index < length && value_array[middle_index].value - value_array[i].value < radius) {
                // Check if the index to the right of the middle index is invalid
                if ((middle_index + 1) >= length || value_array[middle_index + 1].value - value_array[i].value >= radius) {
                    // Found
                    right_index = middle_index;
                } else {
                    // Search right
                    lower_index = middle_index + 1;
                }
            } else {
                // Search left
                higher_index = middle_index - 1;
            }
        }

        total += psi_values[right_index - left_index + 1];
    }

    return total;
}

int partition(Point *point_array, int left_index, int right_index, int coordinate) {
    double pivot_value = point_array[left_index].coords[coordinate];
    int i = left_index + 1;
    int j = right_index;

    while (1) {
        while (point_array[i].coords[coordinate] < pivot_value) {
            i++;
        }

        while (point_array[j].coords[coordinate] > pivot_value) {
            j--;
        }

        if (i >= j) {
            Point temp = point_array[j];
            point_array[j] = point_array[left_index];
            point_array[left_index] = temp;
            return j;
        }

        Point temp = point_array[i];
        point_array[i] = point_array[j];
        point_array[j] = temp;

        i++;
        j--;
    }
}

void partition_median(Point *point_array, int lower_index, int higher_index, int coordinate) {
    int k = (lower_index + higher_index) >> 1;
    while (1) {
        if (lower_index == higher_index) {
            return;
        }

        int pivot_index = partition(point_array, lower_index, higher_index, coordinate);
        if (k == pivot_index) {
            return;
        } else if (k < pivot_index) {
            higher_index = pivot_index - 1;
        } else {
            lower_index = pivot_index + 1;
        }
    }
}


/*
 * axis is the dimension that should be used for this node
 * lower and higher indices are both inclusive and refer to the sorted point array
 * next available index refers to the next index of the kdtree array which can be written to
 *
 * Returns index where it was placed in the kdtree array
 */
int construct_kdtree_helper(KDTree *kdtree_array, Point *point_array, int axis, int lower_index, int higher_index, int *next_available_index) {
    // Select median based on current axis
    int middle_index = (lower_index + higher_index) >> 1;
    double median = point_array[middle_index].coords[axis];

    // Create current node and increment next available tree index
    int selected_index = *next_available_index;
    *next_available_index += 1;
    kdtree_array[selected_index].value = median;
    kdtree_array[selected_index].middle_index = middle_index;

    // Sort again by next axis
    partition_median(point_array, lower_index, middle_index, 1 - axis);
    partition_median(point_array, middle_index + 1, higher_index, 1 - axis);

    // Construct left and right subtrees
    kdtree_array[selected_index].left_index = (lower_index != middle_index) ? construct_kdtree_helper(kdtree_array, point_array, 1 - axis, lower_index, middle_index, next_available_index) : -1;
    kdtree_array[selected_index].right_index = (middle_index + 1 != higher_index) ? construct_kdtree_helper(kdtree_array, point_array, 1 - axis, middle_index + 1, higher_index, next_available_index) : -1;

    return selected_index;
}

void construct_kdtree(KDTree *kdtree_array, Point *point_array, int length) {
    int next_available_index = 0;
    construct_kdtree_helper(kdtree_array, point_array, 0, 0, length - 1, &next_available_index);
}

void update_min_dist_array(double *min_distance_array, int k, double new_min_value) {
    // Search for placement position
    int placement_index = k - 1;
    while (placement_index >= 0 && new_min_value < min_distance_array[placement_index]) {
        placement_index--;
    }
    placement_index++;

    // Update distance rankings
    for (int i = k-1; i > placement_index; i--) {
        min_distance_array[i] = min_distance_array[i-1];
    }
    min_distance_array[placement_index] = new_min_value;
}

void explore_left(KDTree *kdtree_array, Point *point_array, Point point, double *min_distance_array, int k, int axis, int lower_index, int middle_index, int curr_tree_index) {
    // Base case: left subtree is a leaf node
    if (kdtree_array[curr_tree_index].left_index == -1) {
        // Calculate distance between leaf node point and the given point
        double new_dist = max_norm(point_array[lower_index], point);
        if (new_dist < min_distance_array[k-1]) {
            update_min_dist_array(min_distance_array, k, new_dist);
        }
    }

    // Recursive case: left subtree is not a leaf node
    else {
        find_knearest_neighbour_distance_helper(kdtree_array, point_array, point, min_distance_array, k, 1 - axis, lower_index, middle_index, kdtree_array[curr_tree_index].left_index);
    }
}

void explore_right(KDTree *kdtree_array, Point *point_array, Point point, double *min_distance_array, int k, int axis, int middle_index, int higher_index, int curr_tree_index) {
    // Base case: right subtree is a leaf node
    if (kdtree_array[curr_tree_index].right_index == -1) {
        // Calculate distance between leaf node point and the given point
        double new_dist = max_norm(point_array[higher_index], point);
        if (new_dist < min_distance_array[k-1]) {
            update_min_dist_array(min_distance_array, k, new_dist);
        }
    }

    // Recursive case: right subtree is not a leaf node
    else {
        find_knearest_neighbour_distance_helper(kdtree_array, point_array, point, min_distance_array, k, 1 - axis, middle_index + 1, higher_index, kdtree_array[curr_tree_index].right_index);
    }
}

void find_knearest_neighbour_distance_helper(KDTree *kdtree_array, Point *point_array, Point point, double *min_distance_array, int k, int axis, int lower_index, int higher_index, int curr_tree_index) {
    // Use median to decide which subtree to explore first
    int middle_index = kdtree_array[curr_tree_index].middle_index;
    if (point.coords[axis] <= kdtree_array[curr_tree_index].value) {
        // Explore left subtree
        explore_left(kdtree_array, point_array, point, min_distance_array, k, axis, lower_index, middle_index, curr_tree_index);

        // Explore right subtree if needed
        if (kdtree_array[curr_tree_index].value - point.coords[axis] <= min_distance_array[k-1]) {
            explore_right(kdtree_array, point_array, point, min_distance_array, k, axis, middle_index, higher_index, curr_tree_index);
        }
    } else {
        // Explore right subtree
        explore_right(kdtree_array, point_array, point, min_distance_array, k, axis, middle_index, higher_index, curr_tree_index);

        // Explore left subtree if needed
        if (point.coords[axis] - kdtree_array[curr_tree_index].value <= min_distance_array[k-1]) {
            explore_left(kdtree_array, point_array, point, min_distance_array, k, axis, lower_index, middle_index, curr_tree_index);
        }
    }
}

/*
 * Returns the distance between the given point and its kth nearest neighbour in
 * the point array
 */
double find_knearest_neighbour_distance(KDTree *kdtree_array, Point *point_array, int length, Point point, int k) {
    // Initialize min distance array
    double min_distance_array[MAX_K_VALUE];
    for (int i = 0; i < k; i++) {
        min_distance_array[i] = DBL_MAX;
    }

    find_knearest_neighbour_distance_helper(kdtree_array, point_array, point, min_distance_array, k, 0, 0, length - 1, 0);

    return min_distance_array[k-1];
}

/*
 * Normalize the dataset to have zero mean and unit variance.
 */
void normalize_data_set(double *array, int length) {
    // Calculate mean
    double sum = 0;
    for (int i = 0; i < length; i++) {
        sum += array[i];
    }
    double avg = sum / length;

    // Calculate variance
    double sum1 = 0;
    for (int i = 0; i < length; i++) {
        sum1 += pow(array[i] - avg, 2);
    }
    double variance = sum1 / length;

    // Calculate standard deviation
    double std = sqrt(variance);

    // Normalize data set
    for (int i = 0; i < length; i++) {
        array[i] = (array[i] - avg) / std;
    }
}

double _estimate_mi_fast_thread(double *x, double *y, IndexedValue *x_value_array, IndexedValue *y_value_array, double *psi_values, Point *point_array, KDTree *kdtree_array, double *eps, int length, int k) {
    // Construct point array
    for (int i = 0; i < length; i++) {
        point_array[i].coords[0] = x[x_value_array[i].index];
        point_array[i].coords[1] = y[x_value_array[i].index];
    }

    // Consruct KDTree from data points
    construct_kdtree(kdtree_array, point_array, length);

    // Calculate epsilon values
    for (int i = 0; i < length; i++) {
        Point point = { .coords[0] = x[i], .coords[1] = y[i] };
        eps[i] = find_knearest_neighbour_distance(kdtree_array, point_array, length, point, k + 1);
    }

    // Calculate psi average based on x, y, and epsilon
    double psi_avg = (query_ball_point_fast(x_value_array, length, eps, psi_values) + query_ball_point_fast(y_value_array, length, eps, psi_values)) / length;

    return psi_values[length] + psi_values[k] - psi_avg;
}

double _estimate_corr_fast_thread(double *x, double *y, IndexedValue *x_value_array, IndexedValue *y_value_array, double *psi_values, Point *point_array, KDTree *kdtree_array, double *eps, int length, int k) {
    double mi = _estimate_mi_fast_thread(x, y, x_value_array, y_value_array, psi_values, point_array, kdtree_array, eps, length, k);
    return mi <= 0 ? mi : sqrt(1 - exp(-2 * mi));
}

/* --- Public API methods --- */

/*
 * Reads in data from a csv file.
 */
DataShape read_csv(char *filename, double *buffer) {
    // Open file for reading
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Failed to open file\n");
        exit(1);
    }

    // Initialize metadata struct
    DataShape data_shape = { .num_variables = 0, .num_entries_per_variable = 0 };

    // Read file contents
    char *line = NULL;
    size_t len = 0;
    const char delim[2] = ",";
    char *token;
    int i = 0;
    while (getline(&line, &len, file) != -1) {
        token = strtok(line, delim);
        if (token != NULL) {
            sscanf(token, "%lf", buffer + i);
            i++;
        }

        while(token != NULL) {
            token = strtok(NULL, delim);
            if (token != NULL) {
                sscanf(token, "%lf", buffer + i);
                i++;
            }
        }

	// Record number of entries per variable
	if (data_shape.num_entries_per_variable == 0) {
	    data_shape.num_entries_per_variable = i;
	}

	data_shape.num_variables++;
    }

    // Free allocated memory and close the file
    free(line);
    fclose(file);

    return data_shape;
}

/*
 * Saves the given data matrix in a CSV file with the given name.
 */
void save_data_matrix(char *filename, double *data_matrix, int num_variables) {
    // Open file for writing
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Failed to open file\n");
        exit(1);
    }

    // Write contents to the file
    for (int i = 0; i < num_variables; i++) {
        fprintf(file, "%e", data_matrix[i * num_variables]);
        for (int j = 1; j < num_variables; j++) {
            fprintf(file, ",%e", data_matrix[i * num_variables + j]);
        }
        fprintf(file, "\n");
    }
}

int next_var;
pthread_mutex_t next_var_lock;

void *threaded_combine_mi(void *args) {
    // Extract thread arguments
    CombineMiThreadArgs *thread_args = (CombineMiThreadArgs *) args;

    int *valid_values = thread_args->valid_values;
    double *psi_values = thread_args->psi_values;
    double *data = thread_args->data;
    IndexedValue *data_value_arrays = thread_args->data_value_arrays;
    double *data_matrix = thread_args->data_matrix;
    int num_attributes = thread_args->num_attributes;
    int num_variables = thread_args->num_variables;
    int num_entries_per_variable = thread_args->num_entries_per_variable;
    int k = thread_args->k;
    int p = thread_args->p;

    // Allocate heap memory to be used by the thread
    void *total_buffer = malloc(
        num_entries_per_variable * sizeof(Point) +
        num_entries_per_variable * sizeof(KDTree) +
        num_entries_per_variable * sizeof(double)
    );

    if (total_buffer == NULL) {
        fprintf(stderr, "Error allocating memory for GCC computation");
        exit(1);
    }

    Point *point_array = (Point *) total_buffer;
    KDTree *kdtree_array = (KDTree *) (point_array + num_entries_per_variable);
    double *eps = (double *) (kdtree_array + num_entries_per_variable);

    // Represents index into current variable to work on
    int i;

    // Find first variable to work on
    pthread_mutex_lock(&next_var_lock);
    i = next_var;
    if (next_var < num_variables) {
    	next_var++;
    }
    pthread_mutex_unlock(&next_var_lock);

    while (i < num_variables) {
        printf("%d\n", i);
    	data_matrix[i * num_variables + i] = 1;
        for (int j = i + 1; j < num_variables; j++) {
            data_matrix[i * num_variables + j] = 0;
            int num_correlations = 0;

            // Compare all combinations of attributes between variables i and j
            for (int ai = 0; ai < num_attributes; ai++) {
                for (int aj = 0; aj < num_attributes; aj++) {
                    // If attributes are valid for both variables perform correlation calculation
                    if (valid_values[ai * num_variables + i] && valid_values[aj * num_variables + j]) {
                        // Create pointers to the current attribute buffers
                        double *attribute_data_i = data + (ai * num_variables * num_entries_per_variable);
                        double *attribute_data_j = data + (aj * num_variables * num_entries_per_variable);
                        IndexedValue *attribute_data_i_value_array = data_value_arrays + (ai * num_variables * num_entries_per_variable);
                        IndexedValue *attribute_data_j_value_array = data_value_arrays + (aj * num_variables * num_entries_per_variable);

                        // Create pointers to the current variable buffers
                        double *variable_data_i = attribute_data_i + (i * num_entries_per_variable);
                        double *variable_data_j = attribute_data_j + (j * num_entries_per_variable);
                        IndexedValue *variable_data_i_value_array = attribute_data_i_value_array + (i * num_entries_per_variable);
                        IndexedValue *variable_data_j_value_array = attribute_data_j_value_array + (j * num_entries_per_variable);

                        // Estimate generalized correlation coefficient
                        double corr = _estimate_corr_fast_thread(
                            variable_data_i,
                            variable_data_j,
                            variable_data_i_value_array,
                            variable_data_j_value_array,
                            psi_values,
                            point_array,
                            kdtree_array,
                            eps,
                            num_entries_per_variable,
                            k
                        );

                        // If calculated correlation is negative set it to 0
                        corr = (corr >= 0) ? corr : 0;
                        data_matrix[i * num_variables + j] += pow(corr, p);
                        num_correlations++;
                    }
                }
            }

            // Complete generalized mean calculation based on p
            if (num_correlations > 0) {
                data_matrix[i * num_variables + j] = pow(data_matrix[i * num_variables + j] / num_correlations, (1. / p));
            }

            // Use symmetry to save on computation
            data_matrix[j * num_variables + i] = data_matrix[i * num_variables + j];
        }

        // Find next variable to work on
        pthread_mutex_lock(&next_var_lock);
        i = next_var;
        if (next_var < num_variables) {
            next_var++;
        }
        pthread_mutex_unlock(&next_var_lock);
    }
}

double *threaded_combine_mi_starter(double *data, int num_attributes, int num_variables, int num_entries_per_variable, int k, int p, int num_threads) {
	printf("starting\n");

    // Allocate memory needed for computation
    void *total_buffer = malloc(
        (num_variables * num_variables + MAX_NUM_VALUES) * sizeof(double) +
        (num_attributes * num_variables) * sizeof(int) +
        (num_attributes * num_variables * num_entries_per_variable) * sizeof(IndexedValue)
    );

    if (total_buffer == NULL) {
        fprintf(stderr, "Error allocating memory for GCC computation");
        exit(1);
    }

    double *data_matrix = (double *) total_buffer;
    double *psi_values = data_matrix + (num_variables * num_variables);
    int *valid_values = (int *) (psi_values + MAX_NUM_VALUES);
    IndexedValue *data_value_arrays = (IndexedValue *) (valid_values + (num_attributes * num_variables));

    // Determine valid attributes for each variable
    for (int i = 0; i < num_attributes; i++) {
    	for (int j = 0; j < num_variables; j++) {
            valid_values[(i * num_variables) + j] = 1;
            for (int k = 0; k < num_entries_per_variable; k++) {
                if (isnan(data[(i * num_variables * num_entries_per_variable) + (j * num_entries_per_variable) + k])) {
                    valid_values[(i * num_variables) + j] = 0;
                }
            }
        }
    }

    // Compute psi values
    for (int i = 1; i <= num_entries_per_variable; i++) {
        psi_values[i] = psi(i);
    }

    // Normalize values
    for (int i = 0; i < num_attributes; i++) {
        // Create pointers into the current attribute data buffer
        int attribute_data_offset = i * num_variables * num_entries_per_variable;
        double *attribute_data = data + attribute_data_offset;

        for (int j = 0; j < num_variables; j++) {
            // Create pointers into the current variable data buffer
            int variable_data_offset = j * num_entries_per_variable;
            double *variable_data = attribute_data + variable_data_offset;

            // Normalize data set
            normalize_data_set(variable_data, num_entries_per_variable);
        }
    }

    // Sorted arrays compute once
    for (int i = 0; i < num_attributes; i++) {
        // Calculate offset into the current attribute data buffer
        int attribute_data_offset = i * num_variables * num_entries_per_variable;

        // Create pointers to the current data buffer to simplify and improve readability
        double *attribute_data = data + attribute_data_offset;
        IndexedValue *attribute_data_value_array = data_value_arrays + attribute_data_offset;

        // Initialize data
        for (int j = 0; j < num_variables; j++) {
            // Calculate offset into the current variable data buffer
            int variable_data_offset = j * num_entries_per_variable;

            // Create pointers to the current variable buffer to simplify and improve readability
            double *variable_data = attribute_data + variable_data_offset;
            IndexedValue *variable_data_value_array = attribute_data_value_array + variable_data_offset;

            // Initialize indexed data
            for (int k = 0; k < num_entries_per_variable; k++) {
                variable_data_value_array[k].value = variable_data[k];
                variable_data_value_array[k].index = k;
            }

            // Sort initialized data arrays in ascending order by value
            qsort(variable_data_value_array, num_entries_per_variable, sizeof(IndexedValue), compare_value);
        }
    }

    // Initialize next var lock and value
    pthread_mutex_init(&next_var_lock, NULL);
    next_var = 0;

    // Initialize thread arguments
    CombineMiThreadArgs thread_args = {
        .valid_values = valid_values,
        .psi_values = psi_values,
        .data = data,
        .data_value_arrays = data_value_arrays,
        .data_matrix = data_matrix,
        .num_attributes = num_attributes,
        .num_variables = num_variables,
        .num_entries_per_variable = num_entries_per_variable,
        .k = k,
        .p = p
    };

    printf("spinning\n");

    // Spin up threads
    pthread_t tid[MAX_VARIABLES];
    for (int i = 0; i < num_threads; i++) {
    	pthread_create(&(tid[i]), NULL, &threaded_combine_mi, (void *) &thread_args);
    }

    // Wait for all threads to finish
    for (int i = 0; i < num_threads; i++) {
    	pthread_join(tid[i], NULL);
    }

    // Destroy next var lock
    pthread_mutex_destroy(&next_var_lock);

    return data_matrix;
}

/* -- Main -- */

/*
 * First argument should be a comma separated sequence of paths to data files which are csv's of shape
 * (number of variables x number of entries per variable)
 *
 * Second argument should be the output file
 */
int main(int argc, char *argv[]) {
    // Checks for correct number of arguments
    if (argc != 5) {
        printf("Expected 5 arguments\n");
        exit(1);
    }

    // Count number of files given
    int num_files = 1;
    int file_str_len = strlen(argv[1]);
    for (int i = 0; i < file_str_len; i++) {
        if (argv[1][i] == ',') {
            num_files++;
        }
    }

    // Allocate buffer for reading input data
    void *input_buffer = malloc((MAX_VARIABLES * MAX_NUM_VALUES * num_files) * sizeof(double) + (file_str_len + 1) * sizeof(char));
    if (input_buffer == NULL) {
        fprintf(stderr, "Error allocating memory for reading csv");
        exit(1);
    }

    double *data_buffer = (double *) input_buffer;
    char *file_name_buffer = (char *) (data_buffer + (MAX_VARIABLES * MAX_NUM_VALUES * num_files));

    // Read data files
    DataShape data_shape = { .num_variables=0, .num_entries_per_variable=0 };
    int write_index = 0;
    int file_counter = 0;
    for (int i = 0; i < file_str_len + 1; i++) {
        if (argv[1][i] == ',' || argv[1][i] == '\0') {
            // Terminate current file name, read contents, and reset write counter
            file_name_buffer[write_index] = '\0';
            data_shape = read_csv(file_name_buffer, data_buffer + (data_shape.num_variables * data_shape.num_entries_per_variable * file_counter));
            write_index = 0;
            file_counter++;
        } else {
            // Write byte to current file name
            file_name_buffer[write_index] = argv[1][i];
            write_index++;
        }
    }

    printf("Loaded\n");

    // Extract k and p values
    int k = atoi(argv[3]);
    int p = atoi(argv[4]);

    // Estimate MI
    double *data_matrix = threaded_combine_mi_starter(data_buffer, num_files, data_shape.num_variables, data_shape.num_entries_per_variable, k, p, 80);

    // Write data to file
    save_data_matrix(argv[2], data_matrix, data_shape.num_variables);

    // Free allocated memory
    free(input_buffer);
    free(data_matrix);

    return 0;
}
