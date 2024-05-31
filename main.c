#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "rref.h"

typedef struct {
    char *items;
    uint32_t count;
} generic_list_t;

typedef struct {
    double *coefs;
    uint32_t count;
} polynomial_t;

typedef struct {
    polynomial_t *factors;
    uint32_t count;
} factored_t;

typedef struct {
    factored_t *factoreds;
    uint32_t count;
} factored_list_t;

typedef struct {
    polynomial_t *polynomials;
    uint32_t count;
} polynomial_list_t;

void abort_(char *msg) {
    fprintf(stderr, "%s\n", msg);
    exit(1);
}

uint32_t polynomial_coef_count(polynomial_t *p) {
    for (uint32_t i = p->count-1; i >= 0; i--) {
        if (!is_zero(p->coefs[i])) return i+1;
    }
    return 0;
}

const char *superscript_digits[] = {"⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"};

void print_exponent_num(int num) {
    int len = snprintf(NULL, 0, "%d", num);
    char str[len+1];
    snprintf(str, len + 1, "%d", num);
    for (int i = 0; i < len; i++) {
        printf("%s", superscript_digits[str[i] - 0x30]);
    }
}

void print_polynomial(polynomial_t *p) {
    uint32_t coef_count = polynomial_coef_count(p);
    int is_first = 1;
    for (int32_t i = coef_count - 1; i >= 0; i--) {
        double coef = p->coefs[i];
        if (is_first) {
            is_first = 0;
            if (is_double_eq(coef, -1)) {
                printf("-");
                coef = 1;
            }
        } else {
            if (coef < 0) {
                coef *= -1;
                printf(" - ");
            } else {
                printf(" + ");
            }
        }
        if (!is_double_eq(coef, 1) || i == 0) printf("%f", coef);
        if (i > 0) {
            printf("x");
            if (i > 1) {
                print_exponent_num(i);
            }
        }
    }
}

void print_factored(factored_t *f) {
    for (int i = 0; i < f->count; i++) {
        printf("(");
        print_polynomial(&f->factors[i]);
        printf(")");
    }
}

void print_factored_list(factored_list_t *l) {
    for (int i = 0; i < l->count; i++) {
        print_factored(&l->factoreds[i]);
        printf("\n");
    }
}

void print_polynomial_list(polynomial_list_t *l) {
    for (int i = 0; i < l->count; i++) {
        print_polynomial(&l->polynomials[i]);
        printf("\n");
    }
}

generic_list_t *blank_list() {
    return calloc(1, sizeof(generic_list_t));
}

void free_list(generic_list_t *list) {
    free(list->items);
    free(list);
}

uint32_t list_append(generic_list_t *list, uint32_t cap, size_t item_size, void *item) {
    uint32_t count = list->count++;
    if (count == cap) {
        cap = (cap * 3) / 2 + 1;
        list->items = realloc(list->items, cap * item_size);
    }
    memcpy(list->items + count * item_size, item, item_size);
    return cap;
}

void free_factored(factored_t *factored) {
    for (uint32_t i = 0; i < factored->count; i++) {
        free(factored->factors[i].coefs);
    }
    free_list((generic_list_t*)factored);
}

void free_held_factored(factored_t *factored) {
    for (uint32_t i = 0; i < factored->count; i++) {
        free(factored->factors[i].coefs);
    }
    free(factored->factors);
}

void free_polynomial_list(polynomial_list_t *list) {
    for (uint32_t i = 0; i < list->count; i++) {
        free(list->polynomials[i].coefs);
    }
    free_list((generic_list_t*)list);
}

void free_factored_list_shallow(factored_list_t *list) {
    for (uint32_t i = 0; i < list->count; i++) {
        free(list->factoreds[i].factors);
    }
    free_list((generic_list_t*)list);
}

polynomial_t *make_polynomial(double coefs[], int32_t count) {
    polynomial_t *result = malloc(sizeof(polynomial_t));
    double *coefs_mem = malloc(sizeof(double) * count);
    memcpy(coefs_mem, coefs, sizeof(double) * count);
    result->coefs = coefs_mem;
    result->count = count;
    return result;
}

factored_t *make_factored(polynomial_t *factors[], int32_t count) {
    factored_t *result = malloc(sizeof(polynomial_t));
    polynomial_t *factors_mem = malloc(sizeof(polynomial_t) * count);
    for (uint32_t i = 0; i < count; i++) {
        polynomial_t *factor = factors[i];
        memcpy(factors_mem+i, factor, sizeof(polynomial_t));
        free(factor);
    }
    result->factors = factors_mem;
    result->count = count;
    return result;
}

int polynomial_eq(polynomial_t *a, polynomial_t *b) {
    uint32_t count = a->count;
    if (count != b->count) return 0;
    return memcmp(a->coefs, b->coefs, sizeof(double) * count) == 0;
}

uint32_t _all_factored_combos_append(factored_t *factors, factored_list_t *list, uint32_t cap, uint32_t stack[], uint32_t stack_count) {
    polynomial_t *polynomials = malloc(sizeof(polynomial_t) * stack_count);
    for (uint32_t i = 0; i < stack_count; i++) {
        memcpy(polynomials+i, factors->factors+stack[i], sizeof(polynomial_t));
    }
    factored_t factored = {polynomials, stack_count};
    return list_append((generic_list_t*)list, cap, sizeof(factored_t), &factored);
}

uint32_t _all_factored_combos_recurse(factored_t *factors, factored_list_t *list, uint32_t cap, uint32_t stack[], uint32_t *stack_count_p) {
    uint32_t stack_top = 0;
    uint32_t stack_count = *stack_count_p;
    uint32_t new_stack_count = stack_count + 1;
    *stack_count_p = new_stack_count;
    if (stack_count > 0) {
        stack_top = stack[stack_count-1] + 1;
    }
    uint32_t factor_count = factors->count;
    for (uint32_t j = stack_top; j < factor_count; j++) {
        stack[stack_count] = j;
        cap = _all_factored_combos_append(
            factors, list, cap, stack, new_stack_count);
        if (new_stack_count < factor_count - 1) {
            cap = _all_factored_combos_recurse(
                factors, list, cap, stack, stack_count_p);
        }
    }
    (*stack_count_p)--;
    return cap;
}

factored_list_t *all_factored_combos(factored_t *factors) {
    factored_list_t *list = (factored_list_t*)blank_list();
    uint32_t cap = 0;
    uint32_t stack[factors->count-1];
    uint32_t stack_count = 0;
    _all_factored_combos_recurse(factors, list, cap, stack, &stack_count);
    return list;
}

polynomial_t *multiply_polynomials(polynomial_t *a, polynomial_t *b, polynomial_t *result) {
    if (result == NULL) result = malloc(sizeof(polynomial_t));
    uint32_t size = polynomial_coef_count(a) + polynomial_coef_count(b) - 1;
    double *coefs = calloc(1, sizeof(double) * size);
    for (uint32_t i = 0; i < a->count; i++) {
        for (uint32_t j = 0; j < b->count; j++) {
            coefs[i+j] += a->coefs[i] * b->coefs[j];
        }
    }
    result->count = size;
    result->coefs = coefs;
    return result;
}

polynomial_t *add_polynomials(polynomial_t *a, polynomial_t *b, polynomial_t *result) {
    if (result == NULL) result = malloc(sizeof(polynomial_t));
    
    uint32_t ca = polynomial_coef_count(a);
    uint32_t cb = polynomial_coef_count(b);
    
    polynomial_t *first;
    polynomial_t *second;
    uint32_t count;
    uint32_t subcount;
    
    if (ca > cb) {
        count = ca;
        subcount = cb;
        first = a;
        second = b;
    } else {
        count = cb;
        subcount = ca;
        first = b;
        second = a;
    }
    
    double *coefs = malloc(sizeof(double) * count);
    memcpy(coefs+subcount, first->coefs+subcount, sizeof(double) * (count-subcount));
    for (uint32_t i = 0; i < subcount; i++) {
        coefs[i] = first->coefs[i] + second->coefs[i];
    }
    
    result->count = count;
    result->coefs = coefs;
    
    return result;
}

polynomial_t *scale_polynomial(polynomial_t *p, double scale, polynomial_t *result) {
    if (result == NULL) result = malloc(sizeof(polynomial_t));

    uint32_t count = p->count;
    double *coefs = malloc(sizeof(double) * count);
    for (uint32_t i = 0; i < count; i++) {
        coefs[i] = p->coefs[i] * scale;
    }

    result->count = count;
    result->coefs = coefs;

    return result;
}

void expand_factored(factored_t *f, polynomial_t *result) {
    if (f->count == 0) {
        abort_("Cannot expand factored_t with no factors");
    } else if (f->count == 1) {
        uint32_t coef_count = f->factors->count;
        result->count = coef_count;
        size_t coef_mem = sizeof(double) * coef_count;
        double *coefs = malloc(coef_mem);
        memcpy(coefs, f->factors->coefs, coef_mem);
        result->coefs = coefs;
    } else {
        *result = *f->factors;
        for (uint32_t i = 1; i < f->count; i++) {
            double *coefs_mem = result->coefs;
            multiply_polynomials(result, &f->factors[i], result);
            if (i >= 2) free(coefs_mem);
        }
    }
}

polynomial_list_t *expand_factored_list(factored_list_t *list) {
    polynomial_list_t *result = malloc(sizeof(factored_list_t));
    polynomial_t *polynomials = malloc(sizeof(polynomial_t) * list->count);
    result->polynomials = polynomials;
    result->count = list->count;
    for (uint32_t i = 0; i < list->count; i++) {
        expand_factored(&list->factoreds[i], &polynomials[i]);
    }
    return result;
}

void dedup_factored_polynomial_lists(factored_list_t *fi, polynomial_list_t *pi, factored_list_t **fn_p, polynomial_list_t **pn_p) {
    uint32_t count = fi->count;
    if (count != pi->count) {
        abort_("Cannot perform a combined dedup on lists of different lengths");
    }
    
    factored_list_t *fn = malloc(sizeof(factored_list_t));
    factored_t *fs = malloc(sizeof(factored_t) * count);
    fn->factoreds = fs;
    
    polynomial_list_t *pn = malloc(sizeof(polynomial_list_t));
    polynomial_t *ps = malloc(sizeof(polynomial_t) * count);
    pn->polynomials = ps;
    
    uint32_t idx = 0;
    for (uint32_t i = 0; i < count; i++) {
        factored_t *factored = &fi->factoreds[i];
        polynomial_t *polynomial = &pi->polynomials[i];
        for (int32_t j = i-1; j >= 0; j--) {
            polynomial_t *polynomial2 = &pi->polynomials[j];
            if (polynomial_eq(polynomial, polynomial2)) {
                goto is_duplicate;
            }
        }

        fs[idx] = *factored;
        ps[idx] = *polynomial;
        idx++;
        continue;
        
    is_duplicate:
        free(factored->factors);
        free(polynomial->coefs);
        continue;
    }
    
    fn->count = idx;
    pn->count = idx;

    free_list((generic_list_t*)*fn_p);
    free_list((generic_list_t*)*pn_p);
    
    *fn_p = fn;
    *pn_p = pn;
}

factored_list_t *factored_over_factored_list(factored_t *factors, factored_list_t *list) {
    factored_list_t *result = malloc(sizeof(factored_list_t));
    factored_t *factoreds = malloc(sizeof(factored_t) * list->count);
    result->factoreds = factoreds;
    result->count = list->count;
    
    for (uint32_t i = 0; i < list->count; i++) {
        factored_t old = list->factoreds[i];
        factored_t *new_ = factoreds + i;
        uint32_t count = factors->count - old.count;
        new_->count = count;
        polynomial_t *new_factors = malloc(sizeof(polynomial_t) * count);
        new_->factors = new_factors;
        
        uint32_t idx = 0;

        char used_factors[factors->count];
        memset(used_factors, 0, factors->count);
        for (uint32_t j = 0; j < factors->count; j++) {
            polynomial_t factor = factors->factors[j];
            for (uint32_t k = 0; k < old.count; k++) {
                if (used_factors[k]) continue;
                polynomial_t old_factor = old.factors[k];
                if (polynomial_eq(&old_factor, &factor)) {
                    used_factors[k] = 1;
                    goto not_included;
                }
            }

            new_factors[idx++] = factor;

        not_included:;
        }
    }

    return result;
}

void make_matrix(double matrix[], uint32_t matrix_width, uint32_t matrix_height, polynomial_list_t *polynomial_list, polynomial_t *numerator) {
    for (uint32_t x = 0; x < matrix_width - 1; x++) {
        double *cell_p = matrix + x;
        polynomial_t polynomial = polynomial_list->polynomials[x];
        uint32_t y = 0;
        for (; y < polynomial.count; y++) {
            *cell_p = polynomial.coefs[y];
            cell_p += matrix_width;
        }
        for (; y < matrix_height; y++) {
            *cell_p = 0;
            cell_p += matrix_width;
        }
    }
    double *right_p = matrix + (matrix_width - 1);
    for (uint32_t y = 0; y < matrix_height; y++) {
        *right_p = numerator->coefs[y];
        right_p += matrix_width;
    }
}

int extract_leading_values(double matrix[], uint32_t matrix_width, uint32_t matrix_height, double multiples[], uint32_t polynomial_count) {
    memset(multiples, 0, sizeof(double) * polynomial_count);
    uint32_t x = 0;
    double *row = matrix;
    for (uint32_t y = 0; y < matrix_height; y++) {
        double row_end = *(row + (matrix_width - 1));
        for (; x < polynomial_count; x++) {
            if (is_double_eq(*(row+x), 1)) {
                multiples[x] = row_end;
                goto found_leading;
            }
        }
        if (!is_zero(row_end)) {
            return 1;
        }
    found_leading:
        row += matrix_width;
    }
    return 0;
}

void scale_polynomials(polynomial_list_t *polynomials, double multiples[]) {
    for (uint32_t i = 0; i < polynomials->count; i++) {
        polynomial_t *polynomial = &polynomials->polynomials[i];
        double *old_coefs = polynomial->coefs;
        scale_polynomial(polynomial, multiples[i], polynomial);
        free(old_coefs);
    }
}

void filter_zero_multiple_polynomial_list(polynomial_list_t *polynomials, double multiples[]) {
    uint32_t idx = 0;
    for (uint32_t i = 0; i < polynomials->count; i++) {
        polynomial_t polynomial = polynomials->polynomials[i];
        double multiple = multiples[i];
        if (is_zero(multiple)) {
            free(polynomial.coefs);
        } else {
            multiples[idx] = multiple;
            polynomials->polynomials[idx] = polynomial; 
            idx++;
        }
    }
    polynomials->count = idx;
}

void print_decomposed_result(polynomial_list_t *polynomials, double multiples[]) {
    for (uint32_t i = 0; i < polynomials->count; i++) {
        double multiple = multiples[i];
        if (i > 0) {
            if (multiple < 0) {
                printf(" - ");
                multiple *= -1;
            } else {
                printf(" + ");
            }
        }
        printf("%f/(", multiple);
        print_polynomial(&polynomials->polynomials[i]);
        printf(")");
    }
}

int main() {
    polynomial_t *numerator = make_polynomial(
        (double[]) {-44, -8, -26, -2, -4}, 5);
    factored_t *denominator = make_factored((polynomial_t*[]) {
        make_polynomial((double[]) {1, 1}, 2),
        make_polynomial((double[]) {3, 0, 1}, 3),
        make_polynomial((double[]) {3, 0, 1}, 3)
    }, 3);

    printf("(");
    print_polynomial(numerator);
    printf(")/");
    print_factored(denominator);
    printf("\n");

    factored_list_t *factors_list = all_factored_combos(denominator);
    polynomial_list_t *polynomial_list = expand_factored_list(factors_list);
    dedup_factored_polynomial_lists(
        factors_list, polynomial_list, &factors_list, &polynomial_list);

    uint32_t matrix_width = polynomial_list->count + 1;
    uint32_t matrix_height = numerator->count;

    double matrix[matrix_width * matrix_height];

    make_matrix(matrix, matrix_width, matrix_height, polynomial_list, numerator);
    
    rref(matrix, matrix_width, matrix_height);

    double multiples[polynomial_list->count];

    int inconsistent = extract_leading_values(matrix, matrix_width, matrix_height, multiples, polynomial_list->count);

    if (inconsistent) {
        printf("Can't find the partial fraction decomposition\n");
        printf("sorry\n");
        return 0;
    }
    
    factored_list_t *inverse_factors = factored_over_factored_list(
        denominator, factors_list);
    polynomial_list_t *inverse_polynomials = expand_factored_list(inverse_factors);

    filter_zero_multiple_polynomial_list(inverse_polynomials, multiples);

    print_decomposed_result(inverse_polynomials, multiples);
    printf("\n");
    
    free_factored_list_shallow(factors_list);
    free_polynomial_list(polynomial_list);
    free_factored_list_shallow(inverse_factors);
    free_polynomial_list(inverse_polynomials);

    free_list((generic_list_t*)numerator);
    free_factored((factored_t*)denominator);

    return 0;
}
