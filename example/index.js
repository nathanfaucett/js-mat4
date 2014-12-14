global.mat4 = require("../src/index");


var a = mat4.makeRotationX(mat4.create(), Math.PI * 0.25),
    b = mat4.makeRotationZ(mat4.create(), Math.PI * -0.5),
    c = mat4.create();

mat4.mul(c, a, b);

var p = [5, 2, -1],
    s = [2, 0.5, 1],
    r = [0, 0, 0, 1];

mat4.compose(c, p, s, r);
mat4.decompose(c, p, s, r);

console.log(mat4.str(c), p, s, r);
