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
