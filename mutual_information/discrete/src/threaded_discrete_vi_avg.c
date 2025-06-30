#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <pthread.h>

/* Struct definitions */
typedef struct combine_vi_thread_args {
    double *data;
    double *self_entropy_table;
    double *data_matrix;
    int num_attributes;
    int num_rows;
    int num_cols;
    double bin_size;
    int num_bins;
} CombineViThreadArgs;

typedef struct csv {
    int num_rows;
    int num_cols;
    double *buffer;
} CSV;

/*
 * Computes the shannon entropy of x using discretization.
 */
double get_shannon_entropy(double *x, int len, int num_bins, double bin_size, int *histogram_buffer) {
    // Count occurances of each bin value
    memset(histogram_buffer, 0, num_bins * sizeof(int));
    for (int i = 0; i < len; i++) {
        int index = (int) floor(x[i] / bin_size);
        histogram_buffer[index] += 1;
    }

    // Compute shannon entropy
    double entropy = 0;
    for (int i = 0; i < num_bins; i++) {
        if (histogram_buffer[i] != 0) {
            entropy -= histogram_buffer[i] * log(histogram_buffer[i]);
        }
    }

    return (entropy / len) + log(len);
}

/*
 * Computes the joint entropy of x and y using discretization.
 */
double get_joint_entropy(double *x, double *y, int len, int num_bins, double bin_size, int *joint_histogram_buffer, int joint_histogram_buffer_len) {
    // Count occurrances of each bin value
    memset(joint_histogram_buffer, 0, joint_histogram_buffer_len * sizeof(int));
    for (int i = 0; i < len; i++) {
        int x_index = (int) floor(x[i] / bin_size);
        int y_index = (int) floor(y[i] / bin_size);
        if (x_index <= y_index) {
            int index = num_bins * x_index - ((x_index * (x_index + 1)) >> 1) + y_index;
            joint_histogram_buffer[index] += 1;
        } else {
            int index = num_bins * y_index - ((y_index * (y_index + 1)) >> 1) + x_index;
            joint_histogram_buffer[index] += 1;
        }
    }

    // Compute joint entropy
    double entropy = 0;
    for (int i = 0; i < joint_histogram_buffer_len; i++) {
        if (joint_histogram_buffer[i] != 0) {
            entropy -= joint_histogram_buffer[i] * log(joint_histogram_buffer[i]);
        }
    }

    return (entropy / len) + log(len);
}

int next_var;
pthread_mutex_t next_var_lock;

/*
 * This function is used by each thread for computing variation of information
 * values and writing it to a shared data matrix.
 */
void *threaded_combine_vi(void *args) {
    // Extract thread arguments
    CombineViThreadArgs *thread_args = (CombineViThreadArgs *) args;
    double *data = thread_args->data;
    double *self_entropy_table = thread_args->self_entropy_table;
    double *data_matrix = thread_args->data_matrix;
    int num_attributes = thread_args->num_attributes;
    int num_rows = thread_args->num_rows;
    int num_cols = thread_args->num_cols;
    double bin_size = thread_args->bin_size;
    int num_bins = thread_args->num_bins;

    // Allocate heap memory to be used by the thread
    int joint_histogram_buffer_len = (num_bins * (num_bins + 1)) >> 1;
    int *joint_histogram_buffer = malloc(joint_histogram_buffer_len * sizeof(int));
    if (joint_histogram_buffer == NULL) {
        fprintf(stderr, "Error allocating memory for VI computation");
        exit(1);
    }

    // Index into the current variable to work on
    int i;

    // Find first variable to work on
    pthread_mutex_lock(&next_var_lock);
    i = next_var;
    if (next_var < num_rows) {
        next_var++;
    }
    pthread_mutex_unlock(&next_var_lock);

    // Compute VI values
    while (i < num_rows) {
        printf("%d\n", i);
        for (int j = i; j < num_rows; j++) {
            // Initialize information distance to infinity
            data_matrix[i * num_rows + j] = INFINITY;
            int num_distances = 0;

            // Compare all combinations of attributes between variables i and j
            for (int ai = 0; ai < num_attributes; ai++) {
                for (int aj = 0; aj < num_attributes; aj++) {
                    // Get self entropy of each variable
                    double h_x = self_entropy_table[i * num_attributes + ai];
                    double h_y = self_entropy_table[j * num_attributes + aj];

                    // If attributes are valid for both variables perform correlation calculation
                    if (h_x != -1 && h_y != -1) {
                        // Calculate generalized coefficient
                        double *x = data + (ai * num_rows * num_cols) + (i * num_cols);
                        double *y = data + (aj * num_rows * num_cols) + (j * num_cols);
                        double h_x_y = get_joint_entropy(x, y, num_cols, num_bins, bin_size, joint_histogram_buffer, joint_histogram_buffer_len);

                        double vi = 2 * h_x_y - h_x - h_y;

                        // Sum up information distance between the variables
                        data_matrix[i * num_rows + j] += vi;
                        num_distances += 1;
                    }
                }
            }

            // Calculate average by dividing the total information distance by the number of distances
            if (num_distances > 0) {
                data_matrix[i * num_rows + j] = data_matrix[i * num_rows + j] / num_distances;
            }

            // Use symmetry to save on computation
            data_matrix[j * num_rows + i] = data_matrix[i * num_rows + j];
        }

        // Find next variable to work on
        pthread_mutex_lock(&next_var_lock);
        i = next_var;
        if (next_var < num_rows) {
            next_var++;
        }
        pthread_mutex_unlock(&next_var_lock);
    }
}

/*
 * Performs a multithreaded computation of the variation of information matrix
 * of the given data using discretization.
 */
double *threaded_combine_vi_starter(double *data, int num_attributes, int num_rows, int num_cols, double bin_size, int num_bins, int num_threads) {
    // Allocate memory needed for computation
    void *total_buffer = malloc((num_rows * num_rows + num_rows * num_attributes) * sizeof(double) + (num_bins) * sizeof(int) + (num_threads) * sizeof(pthread_t));
    if (total_buffer == NULL) {
        fprintf(stderr, "Error allocating memory for VI computation");
        exit(1);
    }

    double *data_matrix = (double *) total_buffer;
    double *self_entropy_table = data_matrix + (num_rows * num_rows);
    int *histogram_buffer = (int *) (self_entropy_table + (num_rows * num_attributes));

    pthread_t *tid = (pthread_t*) (histogram_buffer + num_bins);

    // Compute self entropy table
    for (int i = 0; i < num_rows; i++) {
        for (int ai = 0; ai < num_attributes; ai++) {
            // Compute current index into the self entropy table
            int table_index = i * num_attributes + ai;

            // Check if the variable contains valid data on this attribute
            self_entropy_table[table_index] = 0;
            for (int j = 0; j < num_cols; j++) {
                if (isnan(data[(ai * num_rows * num_cols) + (i * num_cols) + j])) {
                    self_entropy_table[table_index] = -1;
                }
            }

            // Compute self entropy if valid data exists
            if (self_entropy_table[table_index] != -1) {
                double *x = data + (ai * num_rows * num_cols) + (i * num_cols);
                self_entropy_table[table_index] = get_shannon_entropy(x, num_cols, num_bins, bin_size, histogram_buffer);
            }
        }
    }

    // Initialize locks and next variable value
    pthread_mutex_init(&next_var_lock, NULL);
    next_var = 0;

    // Initialize thread arguments
    CombineViThreadArgs thread_args = {
        .data = data,
        .self_entropy_table = self_entropy_table,
        .data_matrix = data_matrix,
        .num_attributes = num_attributes,
        .num_rows = num_rows,
        .num_cols = num_cols,
        .bin_size = bin_size,
        .num_bins = num_bins
    };

    // Spin up threads
    for (int i = 0; i < num_threads; i++) {
        pthread_create(&(tid[i]), NULL, &threaded_combine_vi, (void *) &thread_args);
    }

    // Wait for all threads to finish
    for (int i = 0; i < num_threads; i++) {
        pthread_join(tid[i], NULL);
    }

    // Destroy next variable lock
    pthread_mutex_destroy(&next_var_lock);

    return data_matrix;
}

/*
 * Parses the string as a positive integer.
 */
int parse_positive_int(char *str, char *value_name) {
    char *endptr;
    long value = strtol(str, &endptr, 10);
    if (endptr == str || *endptr != '\0') {
        fprintf(stderr, "Error parsing %s as an integer\n", value_name);
        exit(1);
    } else if (value <= 0) {
        fprintf(stderr, "%s must be positive");
        exit(1);
    }

    return (int) value;
}

/*
 * Parses thegiven string as a positive floating point number.
 */
double parse_positive_double(char *str, char *value_name) {
    char *endptr;
    double value = strtod(str, &endptr);
    if (value <= 0) {
        fprintf(stderr, "Error parsing %s as a positive floating point value", value_name);
        exit(1);
    }

    return value;
}

/*
 * Reads in data from a csv file.
 */
CSV read_csv(char *filename, double *buffer) {
    // Open file for reading
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Failed to open file\n");
        exit(1);
    }

    // Initialize metadata struct
    CSV csv = { .num_rows = 0, .num_cols = 0, .buffer = buffer };

    // Read file contents
    char *line = NULL;
    size_t len = 0;
    const char delim[2] = ",";
    char *token;
    int i = 0;
    while (getline(&line, &len, file) != -1) {
        token = strtok(line, delim);
        if (token != NULL) {
            sscanf(token, "%lf", csv.buffer + i);
            i++;
        }

        while (token != NULL) {
            token = strtok(NULL, delim);
            if (token != NULL) {
                sscanf(token, "%lf", csv.buffer + i);
                i++;
            }
        }

        // Record number of columns
        if (csv.num_cols == 0) {
            csv.num_cols = i;
        }

        csv.num_rows += 1;
    }

    // Free allocated memory and close the file
    free(line);
    fclose(file);

    return csv;
}

/*
 * Saves the csv to a file with the given name
 */
void save_csv(char *filename, CSV csv) {
    // Open file for writing
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Failed to open file\n");
        exit(1);
    }

    // Write contents to the file
    for (int i = 0; i < csv.num_rows; i++) {
        fprintf(file, "%e", csv.buffer[i * csv.num_rows]);
        for (int j = 1; j < csv.num_cols; j++) {
            fprintf(file, ",%e", csv.buffer[i * csv.num_rows + j]);
        }
        fprintf(file, "\n");
    }
}

/*
 * Computes a variation of information matrix of the given data using
 * discretization and writes it to a CSV file.
 *
 * The first argument should be a comma separated sequence of paths to CSV files
 * of shape (num_variables, num_entries_per_variable). If a variable does not
 * contain valid values, its row should be filled by 'nan'.
 *
 * The second and third arguments should be integers which are greater than or
 * equal to the number of rows and columns of the csv's respectively.
 *
 * The fourth and fifth arguments should be the bin size and number of bins used
 * for discretizing the data.
 *
 * The sixth argument should be the name of the output file to write the MI
 * matrix to as a CSV.
 *
 * The seventh argument should be the number of threads to use for the
 * computation.
 *
 * Example:
 * ```
 * ./threaded_discrete_vi_avg phi.csv,psi.csv 100 300000 1.0 360 output.csv 4
 * ```
 */
int main(int argc, char *argv[]) {
    // Check for correct number of arguments
    if (argc != 8) {
        printf("Expected 7 arguments besides the program name\n");
        exit(1);
    }

    // Count number of files given
    int num_files = 1;
    int file_str_len = strlen(argv[1]);
    for (int i = 0; i < file_str_len; i++) {
        if (argv[1][i] == ',') {
            num_files += 1;
        }
    }

    // Parse numerical command line arguments
    int num_rows_limit = parse_positive_int(argv[2], "number of rows");
    int num_cols_limit = parse_positive_int(argv[3], "number of columns");

    double bin_size = parse_positive_double(argv[4], "bin size");
    int num_bins = parse_positive_int(argv[5], "number of bins");
    int num_threads = parse_positive_int(argv[7], "number of threads");

    // Allocate buffer for reading input data
    void *input_buffer = malloc((num_rows_limit * num_cols_limit * num_files) * sizeof(double) + (file_str_len + 1) * sizeof(char));
    if (input_buffer == NULL) {
        fprintf(stderr, "Error allocating memory for reading csv");
        exit(1);
    }

    double *data_buffer = (double *) input_buffer;
    char *file_name_buffer = (char *) (data_buffer + (num_rows_limit * num_cols_limit * num_files));

    // Read input csv files
    int write_index = 0;
    int file_counter = 0;
    CSV input_csv = { .num_rows = 0, .num_cols = 0, .buffer = NULL };

    for (int i = 0; i < file_str_len + 1; i++) {
        if (argv[1][i] == ',' || argv[1][i] == '\0') {
            // Terminate current file name, read contents, and reset write counter
            file_name_buffer[write_index] = '\0';
            input_csv = read_csv(file_name_buffer, data_buffer + (input_csv.num_rows * input_csv.num_cols * file_counter));
            printf("file loaded: %s\n", file_name_buffer);
            write_index = 0;
            file_counter++;
        } else {
            // Write byte to current file name
            file_name_buffer[write_index] = argv[1][i];
            write_index++;
        }
    }

    // Estimate VI
    double *data_matrix = threaded_combine_vi_starter(data_buffer, num_files, input_csv.num_rows, input_csv.num_cols, bin_size, (int) num_bins, num_threads);

    // Write data to file
    CSV output_csv = { .num_rows = input_csv.num_rows, .num_cols = input_csv.num_rows, .buffer = data_matrix };
    save_csv(argv[6], output_csv);

    // Free allocated memory
    free(input_buffer);
    free(data_matrix);

    return 0;
}
