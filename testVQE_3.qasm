OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg meas[3];
ry(-4.2139206589574325) q[0];
ry(-6.199094870317198) q[1];
cz q[0],q[1];
ry(-3.130149572438669) q[2];
cz q[0],q[2];
ry(-0.9571028797639674) q[0];
cz q[1],q[2];
h q[0];
