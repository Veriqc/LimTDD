OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
z q[6];
h q[2];
h q[7];
z q[5];
y q[6];
y q[0];
z q[0];
h q[4];
y q[1];
x q[8];
z q[6];
s q[3];
z q[7];
s q[1];
z q[6];
s q[9];
x q[1];
z q[2];
x q[9];
h q[4];
cx q[6],q[1];
z q[7];
y q[5];
s q[9];
z q[9];
y q[9];
h q[4];
x q[9];
x q[5];
z q[8];
x q[9];
z q[4];
s q[8];
y q[7];
cx q[9],q[3];
cx q[9],q[1];
s q[0];
cx q[0],q[5];
h q[2];
h q[3];
x q[5];
z q[9];
z q[1];
s q[0];
x q[4];
s q[6];
cx q[3],q[8];
y q[5];
s q[8];
h q[6];
z q[0];
y q[9];
x q[6];
h q[7];
x q[6];
z q[7];
s q[5];
cx q[7],q[9];
h q[7];
x q[9];
z q[2];
cx q[6],q[4];
s q[6];
y q[1];
z q[3];
y q[5];
x q[2];
h q[2];
s q[5];
z q[0];
x q[6];
cx q[4],q[0];
s q[5];
cx q[2],q[4];
z q[1];
y q[0];
x q[3];
y q[2];
z q[5];
s q[6];
s q[5];
x q[4];
y q[5];
s q[4];
s q[7];
x q[4];
y q[4];
y q[3];
z q[2];
h q[1];
s q[2];
h q[6];
s q[1];
y q[1];
cx q[4],q[9];
cx q[9],q[6];
cx q[4],q[1];
x q[7];
x q[5];
z q[5];
x q[7];
h q[0];
s q[4];
x q[5];
cx q[7],q[9];
z q[7];
cx q[0],q[6];
y q[8];
y q[8];
cx q[8],q[1];
y q[4];
cx q[6],q[8];
cx q[7],q[6];
cx q[3],q[7];
y q[6];
h q[2];
s q[5];
x q[7];
z q[0];
y q[8];
cx q[6],q[3];
cx q[2],q[3];
x q[9];
h q[7];
h q[4];
z q[7];
h q[5];
x q[7];
h q[8];
y q[9];
s q[7];
x q[7];
s q[9];
s q[8];
z q[9];
x q[8];
h q[7];
h q[3];
z q[8];
h q[4];
x q[3];
h q[8];
h q[5];
x q[2];
z q[6];
cx q[5],q[0];
h q[9];
y q[3];
z q[0];
z q[9];
cx q[0],q[8];
s q[0];
h q[4];
cx q[4],q[5];
x q[2];
z q[8];
z q[5];
cx q[1],q[9];
s q[8];
cx q[7],q[1];
z q[4];
y q[7];
s q[8];
cx q[5],q[4];
s q[8];
z q[5];
cx q[3],q[2];
z q[9];
x q[1];
s q[1];
x q[5];
h q[0];
s q[7];
x q[8];
cx q[6],q[0];
y q[2];
z q[7];
x q[9];
z q[0];
y q[4];
s q[7];
z q[3];
s q[8];
x q[6];
s q[1];
s q[4];
h q[2];
h q[9];
h q[9];
s q[7];
h q[0];
cx q[8],q[2];
cx q[1],q[8];
h q[1];
s q[7];
y q[9];
y q[7];
z q[1];
h q[8];
h q[2];
