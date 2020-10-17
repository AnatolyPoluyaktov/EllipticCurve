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



}
int main()
{   example();
    //tests ();
    return 0;
}
