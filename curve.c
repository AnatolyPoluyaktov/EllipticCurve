#include <openssl/bn.h>
#include "curve.h"


void InitParam(struct param* param){
    BN_dec2bn(&param->p,p_str);
    BN_dec2bn(&param->a,a_str);
    BN_dec2bn(&param->x_base,x_base_str);
    BN_dec2bn(&param->y_base,y_base_str);
    BN_dec2bn(&param->q, q_str);
    BN_dec2bn(&param->theta, theta_str);

}

void InitJacobiCurve(struct JacobiCurve* curve, struct param* param){
    BN_CTX* tmp = BN_CTX_new ();
    BIGNUM* buf1 = BN_new ();
    BIGNUM* buf2 = BN_new ();
    BN_dec2bn (&curve->Y, "0");
    BN_dec2bn (&curve->X, "0");
    BN_dec2bn (&curve->e, "0");
    BN_dec2bn (&curve->d, "0");
    BN_dec2bn (&curve->Z, "0");
    BN_dec2bn (&curve->p, "0");

    BN_copy(curve->p,param->p);

    //e = -(3 * theta^2 + 4 * a) / 16
    BN_dec2bn (&buf1, "3"); // buf=3
    BN_mod_sqr (curve->e, param->theta, curve->p, tmp); //e = theta ^2
    BN_mod_mul (curve->e, curve->e, buf1, curve->p, tmp);//e = 3 * theta^2
    BN_mod_add (buf2, param->a, param->a, curve->p, tmp); // buf2 = 2*a
    BN_mod_add (buf2, buf2, buf2, curve->p, tmp); //buf2 = 4*a
    BN_mod_add (buf2, curve->e, buf2, curve->p, tmp); //buf2  = 3 * theta^2 + 4 * a
    BN_dec2bn (&buf1, "-1"); //buf1 = -1
    BN_mod_mul (buf2, buf2, buf1, curve->p, tmp); //buf2 = -(3*theta^2+4*a)
    BN_dec2bn (&buf1, "16"); //buf1 = 16
    BN_mod_inverse (buf1, buf1, curve->p, tmp); // buf1 = 1/16
    BN_mod_mul (curve->e, buf1, buf2, curve->p, tmp); //e = -(3 * theta^2 + 4 * a) / 16

    // d = 3 * theta / 4
    BN_mod_add (curve->d, param->theta, param->theta, curve->p, tmp); //d = 2*theta
    BN_mod_add (curve->d, curve->d, param->theta, curve->p, tmp); //d = 3 * theta
    BN_dec2bn (&buf1, "4"); //buf1 = 4
    BN_mod_inverse (buf1, buf1, curve->p, tmp);//buf1 = 1/4
    BN_mod_mul (curve->d, curve->d, buf1, curve->p, tmp);// d = 3 * theta / 4

    // X = 2 * (x_base - theta)
    BN_mod_sub (curve->X, param->x_base, param->theta, curve->p, tmp); // X= x-theta
    BN_mod_add (curve->X, curve->X, curve->X, curve->p, tmp); // X=2*(x - theta)

    // Y = (2 * x_base + theta) * (x_base - theta)^2 - y^2
    BN_mod_add(buf1, param->x_base, param->x_base, curve->p,tmp);// buf1 = 2*x
    BN_mod_add (buf1, buf1,param->theta, curve->p, tmp); //buf1 = 2*x+theta
    BN_mod_sub (buf2, param->x_base, param->theta, curve->p, tmp); //buf2 = x- theta
    BN_mod_sqr (buf2, buf2, curve->p, tmp); // buf2 = (x-theta)^2
    BN_mod_mul (buf2, buf1, buf2, curve->p, tmp); //buf2=(2 * x + theta) * (x - theta)^2
    BN_mod_sqr (curve->Y, param->y_base, curve->p, tmp);// Y= y_base^2
    BN_mod_sub (curve->Y, buf2, curve->Y, curve->p, tmp);//(2 * x + theta) * (x - theta)^2 - y^2

    // Z= y_base
    BN_copy (curve->Z, param->y_base);
    BN_nnmod(curve->Z, curve->Z,curve->p,tmp);

    BN_free(buf1);
    BN_free(buf2);
    BN_CTX_free(tmp);
  }

void InitPoint(struct Point* point, char* X, char* Y, char* Z){
    BN_dec2bn (&point->X, X);
    BN_dec2bn (&point->Y, Y);
    BN_dec2bn (&point->Z, Z);
}
void AdditionPoints(struct Point* P1, struct Point* P2, struct Point* P3, struct JacobiCurve* curve){
    BN_CTX *tmp = BN_CTX_new();

    BIGNUM* T1 = BN_new();
    BN_copy(T1, P1->X);

    BIGNUM* T2 = BN_new();
    BN_copy(T2, P1->Y);

    BIGNUM* T3 = BN_new();
    BN_copy(T3, P1->Z);

    BIGNUM* T4 = BN_new();
    BN_copy(T4, P2->X);

    BIGNUM* T5 = BN_new();
    BN_copy(T5, P2->Y);

    BIGNUM* T6 = BN_new();
    BN_copy(T6, P2->Z);

    BIGNUM* T7 = BN_new();
    BIGNUM* T8 = BN_new();

    //алгоритм сложения

    BN_mod_mul (T7, T1, T3, curve->p, tmp);// T7 = X1 * Z1
    BN_mod_add (T7, T2, T7, curve->p, tmp); // T7 = X1 * Z1 + Y1
    BN_mod_mul (T8, T4, T6, curve->p, tmp); // T8 = X2 * Z2
    BN_mod_add (T8, T5, T8, curve->p, tmp); // T8 = X2 * Z2 + Y2
    BN_mod_mul (T2, T2, T5, curve->p, tmp); // T2 = Y1 * Y2
    BN_mod_mul (T7, T7, T8, curve->p, tmp); // T7 = X3 + Y1 * Y2 + X1 * X2 * Z1 * Z2
    BN_mod_sub (T7, T7, T2, curve->p, tmp); // T7 = X3 + X1 * X2 * Z1 * Z2
    BN_mod_mul (T5, T1, T4, curve->p, tmp); // T5 = X1 * X2
    BN_mod_add (T1, T1, T3, curve->p, tmp); // T1 = X1 + Z1
    BN_mod_mul (T8, T3, T6, curve->p, tmp); // T8 = Z1 * Z2
    BN_mod_add (T4, T4, T6, curve->p, tmp); // T4 = X2 + Z2
    BN_mod_mul (T6, T5, T8, curve->p, tmp); // T6 = X1 * X2 * Z1 * Z2
    BN_mod_sub (T7, T7, T6, curve->p, tmp); // T7 = X3
    BN_mod_mul (T1, T1, T4, curve->p, tmp); // T1 = X1 * Z2 + X2 * Z1 + X1 * X2 + Z1 * Z2
    BN_mod_sub (T1, T1, T5, curve->p, tmp); // T1 = X1 * Z2 + X2 * Z1 + Z1 * Z2
    BN_mod_sub (T1, T1, T8, curve->p, tmp); // T1 = X1 * Z2 + X2 * Z1
    BN_mod_sqr (T3, T1, curve->p, tmp);     // T3 = X1^2 * Z2^2+ X2^2 * Z1^2 + 2 * X1 * X2 * Z1 * Z2
    BN_mod_add (T6, T6, T6, curve->p, tmp); // T6 = 2 * X1 * X2 * Z1 * Z2
    BN_mod_sub (T3, T3, T6, curve->p, tmp); // T3 = X1^2 * Z2^2+ X2^2 * Z1^2
    BN_mod_mul (T4, curve->e, T6, curve->p, tmp);// T4 = 2 * e * X1 * X2 * Z1 * Z2
    BN_mod_mul (T3, T3, T4, curve->p, tmp); // T3 = 2 * e * X1 * X2 * Z1 * Z2 * (X1^2 * Z2^2+ X2^2 * Z1^2)
    BN_mod_mul (T4, curve->d, T6, curve->p, tmp);// T4 = 2 * d * X1 * X2 * Z1 * Z2
    BN_mod_sub (T2, T2, T4, curve->p, tmp); // T2 = Y1 * Y2 - 2 * d * X1 * X2 * Z1 * Z2
    BN_mod_sqr (T4, T8, curve->p, tmp); // T4 = Z1^2 * Z2^2
    BN_mod_sqr (T8, T5, curve->p, tmp);// T8 = X1^2 * X2^2
    BN_mod_mul (T8, curve->e, T8, curve->p, tmp);// T8 = e * X1^2 * X2^2
    BN_mod_add (T5, T4, T8, curve->p, tmp); // T5 = Z1^2 * Z2^2 + e * X1^2 * X2^2
    BN_mod_mul (T2, T2, T5, curve->p, tmp); // T2 = (Z1^2 * Z2^2 + e * X1^2 * X2^2) * (Y1 * Y2 - 2 * d * X1 * X2 * Z1 * Z2)
    BN_mod_add (T2, T2, T3, curve->p, tmp); // T2 = Y3
    BN_mod_sub (T5, T4, T8, curve->p, tmp);// T5 = Z3
    BN_copy (P3->X, T7);  // X3 = T7
    BN_copy (P3->Y, T2); // Y3 = T2
    BN_copy (P3->Z, T5); // Z3 = T5

    BN_CTX_free(tmp);
    BN_free (T1);
    BN_free (T2);
    BN_free (T3);
    BN_free (T4);
    BN_free (T5);
    BN_free (T6);
    BN_free (T7);
    BN_free (T8);

}



void FreeParam(struct param* param){
    BN_free(param->p);
    BN_free(param->a);
    BN_free(param->q);
    BN_free(param->x_base);
    BN_free(param->y_base);

}
void FreeJacobiCurve(struct JacobiCurve* curve){
    BN_free(curve->Y);
    BN_free(curve->X);
    BN_free(curve->e);
    BN_free(curve->d);
    BN_free(curve->Z);
    BN_free(curve->p);
}
