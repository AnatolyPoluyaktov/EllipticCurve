#include <openssl/bn.h>
#include "curve.h"
#include <stdio.h>
#include <string.h>


void example()
{
    struct param prm={NULL, NULL, NULL, NULL, NULL,NULL};
    InitParam(&prm);

    struct JacobiCurve curve={NULL,NULL,NULL,NULL,NULL,NULL};
    InitJacobiCurve(&curve, &prm);
    printf("e=%s\n", BN_bn2dec(curve.e));
    printf("d=%s\n",BN_bn2dec(curve.d));
    printf("X=%s\n",BN_bn2dec(curve.X));
    printf("Y=%s\n",BN_bn2dec(curve.Y));
    printf("Z=%s\n",BN_bn2dec(curve.Z));
    struct Point P = {NULL, NULL, NULL};
    InitPoint(&P,"1","1","1");
    struct Point E = {NULL, NULL, NULL};

    InitPoint(&E,"0","1","1");

    struct Point P3 = {NULL, NULL, NULL};
    InitPoint(&P3,"2","2","2");

    AdditionPoints(&P,&E, &P3,&curve);
    printf("P:\n%s:%s:%s\n",BN_bn2dec(P.X),BN_bn2dec(P.Y),BN_bn2dec(P.Z));
    printf("E:\n%s:%s:%s\n",BN_bn2dec(E.X),BN_bn2dec(E.Y),BN_bn2dec(E.Z));
    printf("P3:\n%s:%s:%s\n",BN_bn2dec(P3.X),BN_bn2dec(P3.Y),BN_bn2dec(P3.Z));







}
int main()
{   example();
    //tests ();
    return 0;
}
