#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//-----some needed const-----
#define SIZE 50
#define SIZEFICT 256
#define MAXCONST -1000000000000
#define MINCONST 1000000000000
#define STEP 0.001
#define APPRVALUE 0.01

//----struct of point, parametric point, spline------
typedef struct{
    double x;
    double y;
} POINT;

typedef struct{
    double t;
    POINT point;
} PARAMETRPOINT;


typedef struct{
    PARAMETRPOINT* mass_points;
    POINT* mass_si;
    int len;
} SPLINE;


//----reading file-----
PARAMETRPOINT* read_file(char name_of_file[], int len){
    FILE *S1;
    int file_iter = 0;
    PARAMETRPOINT* mass_points = malloc(len * sizeof(PARAMETRPOINT));
    S1 = fopen(name_of_file, "r");
    while(fscanf(S1, "%lf%lf%lf", &mass_points[file_iter].t, &mass_points[file_iter].point.x, &mass_points[file_iter].point.y) != EOF){
        file_iter++;
    }
    fclose(S1);
    return mass_points;
}


int count_string_file(char name_of_file[]){
    int ans = 0;
    char fict[SIZEFICT];
    FILE *S1;
    S1 = fopen(name_of_file, "r");
    while(fgets(fict, SIZEFICT-1, S1) != NULL){
        ans++;
    }
    fclose(S1);
    return ans;
}


//---find intersection of array of x betweeen 2 splines-----
POINT intersection(PARAMETRPOINT* first_mass_x, int first_len, PARAMETRPOINT* second_mass_x, int second_len){
    double first_min, first_max, second_min, second_max;
    POINT intersect;
    first_min = MINCONST;
    second_min = MINCONST;
    first_max = MAXCONST;
    first_max = MAXCONST;
    for(int i = 0; i < first_len; i++){
        first_min = (first_min <= first_mass_x[i].point.x) ? first_min : first_mass_x[i].point.x;
        first_max = (first_max >= first_mass_x[i].point.x) ? first_max : first_mass_x[i].point.x;
    }
    for(int i = 0; i < second_len; i++){
        second_min = (second_min <= second_mass_x[i].point.x) ? first_min : second_mass_x[i].point.x;
        second_max = (second_max >= second_mass_x[i].point.x) ? first_max : second_mass_x[i].point.x;
    }
    if(first_min >= second_min){
        if(first_max <= second_max){
            intersect.x = first_min;
            intersect.y = first_max;
            return intersect;
        } else {
            intersect.x = first_min;
            intersect.y = second_max;
            if(first_min > second_max){
                intersect.x = MINCONST;
                intersect.y = MINCONST;
                return intersect;
            }
            return intersect;
        }
    } else {
        if(first_max <= second_max){
            intersect.x = second_min;
            intersect.y = first_max;
            if(second_min > first_max){
                intersect.x = MINCONST;
                intersect.y = MINCONST;
                return intersect;
            }
            return intersect;
        } else {
            intersect.x = second_min;
            intersect.y = second_max;
            return intersect;
        }
    }
}


//---form matrix of si------
double** matrix_of_si(PARAMETRPOINT* mass_points, int len){
    double** diagonal_3_matrix = malloc((len-2)*sizeof(double*));
    for(int i = 0; i < len-2; i++){
        double* row_matrix = calloc(len-2, sizeof(double));
        if(0 == i){
            row_matrix[0] = 2*(mass_points[2].t-mass_points[0].t);
            row_matrix[1] = (mass_points[2].t-mass_points[1].t);
        } else if(len-3 == i){
            row_matrix[i-1] = (mass_points[i+1].t-mass_points[i].t);
            row_matrix[i] = 2*(mass_points[i+2].t-mass_points[i].t);
        } else{
            row_matrix[i-1] = (mass_points[i+1].t-mass_points[i].t);
            row_matrix[i] = 2*(mass_points[i+2].t-mass_points[i].t);
            row_matrix[i+1] = (mass_points[i+2].t-mass_points[i+1].t);
        }
        diagonal_3_matrix[i] = row_matrix;
    }
    return diagonal_3_matrix;
}


//---forming matrix of right equation-----
POINT* matrix_of_p(PARAMETRPOINT* mass_points, int len){
    POINT* p_matrix = malloc((len-2)*sizeof(POINT));
    for(int i = 0; i < len-2; i++){
        p_matrix[i].x = 6*(mass_points[i+2].point.x - mass_points[i+1].point.x)/(mass_points[i+2].t - mass_points[i+1].t) - 6*(mass_points[i+1].point.x-mass_points[i].point.x)/(mass_points[i+1].t - mass_points[i].t);
        p_matrix[i].y = 6*(mass_points[i+2].point.y - mass_points[i+1].point.y)/(mass_points[i+2].t - mass_points[i+1].t) - 6*(mass_points[i+1].point.y-mass_points[i].point.y)/(mass_points[i+1].t - mass_points[i].t);
    }
    return p_matrix;
}

//----solving matrix of si------
POINT* matrix_solution(double** a, POINT* b, int len_a){
    POINT* solution = calloc(len_a+2, sizeof(POINT));
    double* x = calloc(len_a, sizeof(double));
    double* v = calloc(len_a, sizeof(double));
    double* u = calloc(len_a, sizeof(double));

    v[0] = a[0][1] / (-a[0][0]);
    u[0] = (- b[0].x) / (-a[0][0]);

    for(int i = 1; i < len_a-1; i++){
        v[i] = a[i][i + 1] / (-a[i][i] - a[i][i - 1] * v[i - 1]);
        u[i] = (a[i][i - 1] * u[i - 1] - b[i].x) / (-a[i][i] - a[i][i - 1] * v[i - 1]);
    }
    v[len_a - 1] = 0;
    u[len_a - 1] = (a[len_a - 1][len_a - 2] * u[len_a - 2] - b[len_a - 1].x) / (-a[len_a - 1][len_a - 1] - a[len_a - 1][len_a - 2] * v[len_a - 2]);

    x[len_a - 1] = u[len_a - 1];
    for(int i = len_a - 1; i > 0; i--){
        x[i - 1] = v[i - 1] * x[i] + u[i - 1];
    }
    for(int i = 1; i < len_a+1; i++){
        solution[i].x = x[i-1];
        x[i-1] = 0;
    }
    v[0] = a[0][1] / (-a[0][0]);
    u[0] = (- b[0].y) / (-a[0][0]);

    for(int i = 1; i < len_a-1; i++){
        v[i] = a[i][i + 1] / (-a[i][i] - a[i][i - 1] * v[i - 1]);
        u[i] = (a[i][i - 1] * u[i - 1] - b[i].y) / (-a[i][i] - a[i][i - 1] * v[i - 1]);
    }
    v[len_a - 1] = 0;
    u[len_a - 1] = (a[len_a - 1][len_a - 2] * u[len_a - 2] - b[len_a - 1].y) / (-a[len_a - 1][len_a - 1] - a[len_a - 1][len_a - 2] * v[len_a - 2]);

    x[len_a - 1] = u[len_a - 1];
    for(int i = len_a - 1; i > 0; i--){
        x[i - 1] = v[i - 1] * x[i] + u[i - 1];
    }
    for(int i = 1; i < len_a+1; i++){
        solution[i].y = x[i-1];
    }

    return solution;

}

//---forming a spline----
SPLINE form_spline(PARAMETRPOINT* mass_points, int len){
    SPLINE spline;
    double** mass_diagonal = matrix_of_si(mass_points, len);
    POINT* mass_p = matrix_of_p(mass_points, len);
    POINT* matrix_si = matrix_solution(mass_diagonal, mass_p, len-2);
    spline.mass_points = mass_points;
    spline.mass_si = matrix_si;
    spline.len = len;
    return spline;
}

SPLINE create_spline(){
    int len;
    char name[SIZE];
    printf("Print path to file\n");
    scanf("%s", &name);
    len = count_string_file(name);
    PARAMETRPOINT* mass_points = read_file(name, len);
    return form_spline(mass_points, len);
}


//----find a point of spline----
POINT func_spline(SPLINE spline, double t){
    POINT point;
    int i = 0;
    for(int j=0; j < spline.len-1; j++){
        if((t <= spline.mass_points[j+1].t) && (t >= spline.mass_points[j].t)){
            i = j;
            break;
        }
    }
    double w = (t-spline.mass_points[i].t)/(spline.mass_points[i+1].t-spline.mass_points[i].t);
    point.x = (1-w)*spline.mass_points[i].point.x+w*spline.mass_points[i+1].point.x+((-2*w+3*w*w-w*w*w)*spline.mass_si[i].x+(-w+w*w*w)*spline.mass_si[i+1].x)*(spline.mass_points[i+1].t-spline.mass_points[i].t)*(spline.mass_points[i+1].t-spline.mass_points[i].t)/6;
    point.y = (1-w)*spline.mass_points[i].point.y+w*spline.mass_points[i+1].point.y+((-2*w+3*w*w-w*w*w)*spline.mass_si[i].y+(-w+w*w*w)*spline.mass_si[i+1].y)*(spline.mass_points[i+1].t-spline.mass_points[i].t)*(spline.mass_points[i+1].t-spline.mass_points[i].t)/6;
    return point;
}


//---find intersection point----
POINT find_point(POINT intersection_x, SPLINE first_spline, SPLINE second_spline){
    double point_i = 0;
    double point_i_2 = 0;
    POINT point;
    point.x = 0;
    point.y = 0;
    POINT point_2;
    point_2.x = 0;
    point_2.y = 0;
    POINT ans;
    while(point_i<(double)first_spline.len){
        point_i_2 = 0;
        point = func_spline(first_spline, point_i);
        if(point.x<intersection_x.x || point.x>intersection_x.y){
            point_i += STEP;
            continue;
        }
        while(point_i_2<(double)second_spline.len){
            point_2 = func_spline(second_spline, point_i_2);
            if(point_2.x<intersection_x.x || point_2.x>intersection_x.y){
                point_i_2 += STEP;
                continue;
            }
            if((fabs(point_2.x - point.x) < APPRVALUE) && (fabs(point_2.y - point.y) < APPRVALUE)){
                ans.x = point.x;
                ans.y = point.y;
                return ans;
            }
            point_i_2 += STEP;
        }
        point_i += STEP;
    }
    ans.x = NAN;
    ans.y = NAN;
    return ans;
}


//---find distance between splines----
double find_min_distance(SPLINE first_spline, SPLINE second_spline){
    double point_i = 0;
    double point_i_2 = 0;
    POINT point;
    point.x = 0;
    point.y = 0;
    POINT point_2;
    point_2.x = 0;
    point_2.y = 0;
    double min_distance = MINCONST;
    while(point_i<(double)first_spline.len){
        point_i_2 = 0;
        point = func_spline(first_spline, point_i);
        while(point_i_2<(double)second_spline.len){
            point_2 = func_spline(second_spline, point_i_2);
            min_distance = (min_distance <= sqrt(pow(point_2.x-point.x,2)+pow(point_2.y-point.y, 2))) ? min_distance : sqrt(pow(point_2.x-point.x,2)+pow(point_2.y-point.y, 2));
            point_i_2 += STEP;
        }
        point_i += STEP;
    }
    return min_distance;
}


int main() {
    SPLINE first_spline = create_spline();
    SPLINE second_spline = create_spline();

    POINT intersection_x = intersection(first_spline.mass_points, first_spline.len, second_spline.mass_points, second_spline.len);
    POINT point_of_intersection = find_point(intersection_x, first_spline, second_spline);
    printf("\nINTERSECTION OF X: [%lf, %lf]\n", intersection_x.x, intersection_x.y);
    printf("\nPOINT OF INTERSECTION: {%lf; %lf}\n", point_of_intersection.x, point_of_intersection.y);
    printf("\nDISTANCE: %lf\n", find_min_distance(first_spline, second_spline));
    /*double point_i = 0;
    POINT current_point;
    while(point_i < 6){
        current_point = func_spline(first_spline, point_i);
        printf("(%lf; %lf) ", current_point.x, current_point.y);
        point_i = point_i + 0.1;
    }*/
    return 0;
}
