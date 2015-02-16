var mathf = require("mathf"),
    vec3 = require("vec3");


var mat4 = exports;


mat4.ArrayType = typeof(Float32Array) !== "undefined" ? Float32Array : mathf.ArrayType;


mat4.create = function(m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44) {
    var out = new mat4.ArrayType(16);

    out[0] = m11 !== undefined ? m11 : 1;
    out[4] = m12 !== undefined ? m12 : 0;
    out[8] = m13 !== undefined ? m13 : 0;
    out[12] = m14 !== undefined ? m14 : 0;
    out[1] = m21 !== undefined ? m21 : 0;
    out[5] = m22 !== undefined ? m22 : 1;
    out[9] = m23 !== undefined ? m23 : 0;
    out[13] = m24 !== undefined ? m24 : 0;
    out[2] = m31 !== undefined ? m31 : 0;
    out[6] = m32 !== undefined ? m32 : 0;
    out[10] = m33 !== undefined ? m33 : 1;
    out[14] = m34 !== undefined ? m34 : 0;
    out[3] = m41 !== undefined ? m41 : 0;
    out[7] = m42 !== undefined ? m42 : 0;
    out[11] = m43 !== undefined ? m43 : 0;
    out[15] = m44 !== undefined ? m44 : 1;

    return out;
};

mat4.copy = function(out, a) {

    out[0] = a[0];
    out[1] = a[1];
    out[2] = a[2];
    out[3] = a[3];
    out[4] = a[4];
    out[5] = a[5];
    out[6] = a[6];
    out[7] = a[7];
    out[8] = a[8];
    out[9] = a[9];
    out[10] = a[10];
    out[11] = a[11];
    out[12] = a[12];
    out[13] = a[13];
    out[14] = a[14];
    out[15] = a[15];

    return out;
};

mat4.clone = function(a) {
    var out = new mat4.ArrayType(16);

    out[0] = a[0];
    out[1] = a[1];
    out[2] = a[2];
    out[3] = a[3];
    out[4] = a[4];
    out[5] = a[5];
    out[6] = a[6];
    out[7] = a[7];
    out[8] = a[8];
    out[9] = a[9];
    out[10] = a[10];
    out[11] = a[11];
    out[12] = a[12];
    out[13] = a[13];
    out[14] = a[14];
    out[15] = a[15];

    return out;
};

mat4.set = function(out, m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44) {

    out[0] = m11 !== undefined ? m11 : 1;
    out[4] = m12 !== undefined ? m12 : 0;
    out[8] = m13 !== undefined ? m13 : 0;
    out[12] = m14 !== undefined ? m14 : 0;
    out[1] = m21 !== undefined ? m21 : 0;
    out[5] = m22 !== undefined ? m22 : 1;
    out[9] = m23 !== undefined ? m23 : 0;
    out[13] = m24 !== undefined ? m24 : 0;
    out[2] = m31 !== undefined ? m31 : 0;
    out[6] = m32 !== undefined ? m32 : 0;
    out[10] = m33 !== undefined ? m33 : 1;
    out[14] = m34 !== undefined ? m34 : 0;
    out[3] = m41 !== undefined ? m41 : 0;
    out[7] = m42 !== undefined ? m42 : 0;
    out[11] = m43 !== undefined ? m43 : 0;
    out[15] = m44 !== undefined ? m44 : 1;

    return out;
};

mat4.identity = function(out) {

    out[0] = 1;
    out[1] = 0;
    out[2] = 0;
    out[3] = 0;
    out[4] = 0;
    out[5] = 1;
    out[6] = 0;
    out[7] = 0;
    out[8] = 0;
    out[9] = 0;
    out[10] = 1;
    out[11] = 0;
    out[12] = 0;
    out[13] = 0;
    out[14] = 0;
    out[15] = 1;

    return out;
};

mat4.zero = function(out) {

    out[0] = 0;
    out[1] = 0;
    out[2] = 0;
    out[3] = 0;
    out[4] = 0;
    out[5] = 0;
    out[6] = 0;
    out[7] = 0;
    out[8] = 0;
    out[9] = 0;
    out[10] = 0;
    out[11] = 0;
    out[12] = 0;
    out[13] = 0;
    out[14] = 0;
    out[15] = 0;

    return out;
};

mat4.mul = function(out, a, b) {
    var a11 = a[0],
        a12 = a[4],
        a13 = a[8],
        a14 = a[12],
        a21 = a[1],
        a22 = a[5],
        a23 = a[9],
        a24 = a[13],
        a31 = a[2],
        a32 = a[6],
        a33 = a[10],
        a34 = a[14],
        a41 = a[3],
        a42 = a[7],
        a43 = a[11],
        a44 = a[15],

        b11 = b[0],
        b12 = b[4],
        b13 = b[8],
        b14 = b[12],
        b21 = b[1],
        b22 = b[5],
        b23 = b[9],
        b24 = b[13],
        b31 = b[2],
        b32 = b[6],
        b33 = b[10],
        b34 = b[14],
        b41 = b[3],
        b42 = b[7],
        b43 = b[11],
        b44 = b[15];

    out[0] = a11 * b11 + a12 * b21 + a13 * b31 + a14 * b41;
    out[4] = a11 * b12 + a12 * b22 + a13 * b32 + a14 * b42;
    out[8] = a11 * b13 + a12 * b23 + a13 * b33 + a14 * b43;
    out[12] = a11 * b14 + a12 * b24 + a13 * b34 + a14 * b44;

    out[1] = a21 * b11 + a22 * b21 + a23 * b31 + a24 * b41;
    out[5] = a21 * b12 + a22 * b22 + a23 * b32 + a24 * b42;
    out[9] = a21 * b13 + a22 * b23 + a23 * b33 + a24 * b43;
    out[13] = a21 * b14 + a22 * b24 + a23 * b34 + a24 * b44;

    out[2] = a31 * b11 + a32 * b21 + a33 * b31 + a34 * b41;
    out[6] = a31 * b12 + a32 * b22 + a33 * b32 + a34 * b42;
    out[10] = a31 * b13 + a32 * b23 + a33 * b33 + a34 * b43;
    out[14] = a31 * b14 + a32 * b24 + a33 * b34 + a34 * b44;

    out[3] = a41 * b11 + a42 * b21 + a43 * b31 + a44 * b41;
    out[7] = a41 * b12 + a42 * b22 + a43 * b32 + a44 * b42;
    out[11] = a41 * b13 + a42 * b23 + a43 * b33 + a44 * b43;
    out[15] = a41 * b14 + a42 * b24 + a43 * b34 + a44 * b44;

    return out;
};

mat4.smul = function(out, a, s) {

    out[0] = a[0] * s;
    out[1] = a[1] * s;
    out[2] = a[2] * s;
    out[3] = a[3] * s;
    out[4] = a[4] * s;
    out[5] = a[5] * s;
    out[6] = a[6] * s;
    out[7] = a[7] * s;
    out[8] = a[8] * s;
    out[9] = a[9] * s;
    out[10] = a[10] * s;
    out[11] = a[11] * s;
    out[12] = a[12] * s;
    out[13] = a[13] * s;
    out[14] = a[14] * s;
    out[15] = a[15] * s;

    return out;
};

mat4.sdiv = function(out, a, s) {
    s = s !== 0 ? 1 / s : s;

    out[0] = a[0] * s;
    out[1] = a[1] * s;
    out[2] = a[2] * s;
    out[3] = a[3] * s;
    out[4] = a[4] * s;
    out[5] = a[5] * s;
    out[6] = a[6] * s;
    out[7] = a[7] * s;
    out[8] = a[8] * s;
    out[9] = a[9] * s;
    out[10] = a[10] * s;
    out[11] = a[11] * s;
    out[12] = a[12] * s;
    out[13] = a[13] * s;
    out[14] = a[14] * s;
    out[15] = a[15] * s;

    return out;
};

mat4.determinant = function(out) {
    var m11 = out[0],
        m12 = out[4],
        m13 = out[8],
        m14 = out[12],
        m21 = out[1],
        m22 = out[5],
        m23 = out[9],
        m24 = out[13],
        m31 = out[2],
        m32 = out[6],
        m33 = out[10],
        m34 = out[14],
        m41 = out[3],
        m42 = out[7],
        m43 = out[11],
        m44 = out[15];

    return (
        m41 * (m14 * m23 * m32 - m13 * m24 * m32 - m14 * m22 * m33 + m12 * m24 * m33 + m13 * m22 * m34 - m12 * m23 * m34) +
        m42 * (m11 * m23 * m34 - m11 * m24 * m33 + m14 * m21 * m33 - m13 * m21 * m34 + m13 * m24 * m31 - m14 * m23 * m31) +
        m43 * (m11 * m24 * m32 - m11 * m22 * m34 - m14 * m21 * m32 + m12 * m21 * m34 + m14 * m22 * m31 - m12 * m24 * m31) +
        m44 * (-m13 * m22 * m31 - m11 * m23 * m32 + m11 * m22 * m33 + m13 * m21 * m32 - m12 * m21 * m33 + m12 * m23 * m31)
    );
};

mat4.inverse = function(out, a) {
    var m11 = a[0],
        m12 = a[4],
        m13 = a[8],
        m14 = a[12],
        m21 = a[1],
        m22 = a[5],
        m23 = a[9],
        m24 = a[13],
        m31 = a[2],
        m32 = a[6],
        m33 = a[10],
        m34 = a[14],
        m41 = a[3],
        m42 = a[7],
        m43 = a[11],
        m44 = a[15],

        me0 = m23 * m34 * m42 - m24 * m33 * m42 + m24 * m32 * m43 - m22 * m34 * m43 - m23 * m32 * m44 + m22 * m33 * m44,
        me4 = m14 * m33 * m42 - m13 * m34 * m42 - m14 * m32 * m43 + m12 * m34 * m43 + m13 * m32 * m44 - m12 * m33 * m44,
        me8 = m13 * m24 * m42 - m14 * m23 * m42 + m14 * m22 * m43 - m12 * m24 * m43 - m13 * m22 * m44 + m12 * m23 * m44,
        me12 = m14 * m23 * m32 - m13 * m24 * m32 - m14 * m22 * m33 + m12 * m24 * m33 + m13 * m22 * m34 - m12 * m23 * m34,

        det = m11 * me0 + m21 * me4 + m31 * me8 + m41 * me12;

    if (det === 0) {
        return mat4.identity(out);
    }
    det = 1 / det;

    out[0] = me0 * det;
    out[4] = me4 * det;
    out[8] = me8 * det;
    out[12] = me12 * det;
    out[1] = (m24 * m33 * m41 - m23 * m34 * m41 - m24 * m31 * m43 + m21 * m34 * m43 + m23 * m31 * m44 - m21 * m33 * m44) * det;
    out[5] = (m13 * m34 * m41 - m14 * m33 * m41 + m14 * m31 * m43 - m11 * m34 * m43 - m13 * m31 * m44 + m11 * m33 * m44) * det;
    out[9] = (m14 * m23 * m41 - m13 * m24 * m41 - m14 * m21 * m43 + m11 * m24 * m43 + m13 * m21 * m44 - m11 * m23 * m44) * det;
    out[13] = (m13 * m24 * m31 - m14 * m23 * m31 + m14 * m21 * m33 - m11 * m24 * m33 - m13 * m21 * m34 + m11 * m23 * m34) * det;
    out[2] = (m22 * m34 * m41 - m24 * m32 * m41 + m24 * m31 * m42 - m21 * m34 * m42 - m22 * m31 * m44 + m21 * m32 * m44) * det;
    out[6] = (m14 * m32 * m41 - m12 * m34 * m41 - m14 * m31 * m42 + m11 * m34 * m42 + m12 * m31 * m44 - m11 * m32 * m44) * det;
    out[10] = (m12 * m24 * m41 - m14 * m22 * m41 + m14 * m21 * m42 - m11 * m24 * m42 - m12 * m21 * m44 + m11 * m22 * m44) * det;
    out[14] = (m14 * m22 * m31 - m12 * m24 * m31 - m14 * m21 * m32 + m11 * m24 * m32 + m12 * m21 * m34 - m11 * m22 * m34) * det;
    out[3] = (m23 * m32 * m41 - m22 * m33 * m41 - m23 * m31 * m42 + m21 * m33 * m42 + m22 * m31 * m43 - m21 * m32 * m43) * det;
    out[7] = (m12 * m33 * m41 - m13 * m32 * m41 + m13 * m31 * m42 - m11 * m33 * m42 - m12 * m31 * m43 + m11 * m32 * m43) * det;
    out[11] = (m13 * m22 * m41 - m12 * m23 * m41 - m13 * m21 * m42 + m11 * m23 * m42 + m12 * m21 * m43 - m11 * m22 * m43) * det;
    out[15] = (m12 * m23 * m31 - m13 * m22 * m31 + m13 * m21 * m32 - m11 * m23 * m32 - m12 * m21 * m33 + m11 * m22 * m33) * det;

    return out;
};

mat4.transpose = function(out, a) {
    var a01, a02, a03, a12, a13, a23;

    if (out === a) {
        a01 = a[1];
        a02 = a[2];
        a03 = a[3];
        a12 = a[6];
        a13 = a[7];
        a23 = a[11];

        out[1] = a[4];
        out[2] = a[8];
        out[3] = a[12];
        out[4] = a01;
        out[6] = a[9];
        out[7] = a[13];
        out[8] = a02;
        out[9] = a12;
        out[11] = a[14];
        out[12] = a03;
        out[13] = a13;
        out[14] = a23;
    } else {
        out[0] = a[0];
        out[1] = a[4];
        out[2] = a[8];
        out[3] = a[12];
        out[4] = a[1];
        out[5] = a[5];
        out[6] = a[9];
        out[7] = a[13];
        out[8] = a[2];
        out[9] = a[6];
        out[10] = a[10];
        out[11] = a[14];
        out[12] = a[3];
        out[13] = a[7];
        out[14] = a[11];
        out[15] = a[15];
    }

    return out;
};

var lookAt_x = vec3.create(),
    lookAt_y = vec3.create(),
    lookAt_z = vec3.create();

mat4.lookAt = function(out, eye, target, up) {
    var x = lookAt_x,
        y = lookAt_y,
        z = lookAt_z;

    vec3.sub(z, eye, target);
    vec3.normalize(z, z);

    if (vec3.length(z) === 0) {
        z[2] = 1;
    }

    vec3.cross(x, up, z);
    vec3.normalize(x, x);

    if (vec3.length(x) === 0) {
        z[0] += mathf.EPSILON;
        vec3.cross(x, up, z);
        vec3.normalize(x, x);
    }

    vec3.cross(y, z, x);

    out[0] = x[0];
    out[4] = y[0];
    out[8] = z[0];
    out[1] = x[1];
    out[5] = y[1];
    out[9] = z[1];
    out[2] = x[2];
    out[6] = y[2];
    out[10] = z[2];

    return out;
};

mat4.compose = function(out, position, scale, rotation) {
    var x = rotation[0],
        y = rotation[1],
        z = rotation[2],
        w = rotation[3],
        x2 = x + x,
        y2 = y + y,
        z2 = z + z,
        xx = x * x2,
        xy = x * y2,
        xz = x * z2,
        yy = y * y2,
        yz = y * z2,
        zz = z * z2,
        wx = w * x2,
        wy = w * y2,
        wz = w * z2,

        sx = scale[0],
        sy = scale[1],
        sz = scale[2];

    out[0] = (1 - (yy + zz)) * sx;
    out[4] = (xy - wz) * sy;
    out[8] = (xz + wy) * sz;

    out[1] = (xy + wz) * sx;
    out[5] = (1 - (xx + zz)) * sy;
    out[9] = (yz - wx) * sz;

    out[2] = (xz - wy) * sx;
    out[6] = (yz + wx) * sy;
    out[10] = (1 - (xx + yy)) * sz;

    out[3] = 0;
    out[7] = 0;
    out[11] = 0;

    out[12] = position[0];
    out[13] = position[1];
    out[14] = position[2];
    out[15] = 1;

    return out;
};

mat4.decompose = function(out, position, scale, rotation) {
    var m11 = out[0],
        m12 = out[4],
        m13 = out[8],
        m21 = out[1],
        m22 = out[5],
        m23 = out[9],
        m31 = out[2],
        m32 = out[6],
        m33 = out[10],
        x = 0,
        y = 0,
        z = 0,
        w = 1,

        sx = vec3.lengthValues(m11, m21, m31),
        sy = vec3.lengthValues(m12, m22, m32),
        sz = vec3.lengthValues(m13, m23, m33),

        invSx = 1 / sx,
        invSy = 1 / sy,
        invSz = 1 / sz,

        s, trace;

    scale[0] = sx;
    scale[1] = sy;
    scale[2] = sz;

    position[0] = out[12];
    position[1] = out[13];
    position[2] = out[14];

    m11 *= invSx;
    m12 *= invSy;
    m13 *= invSz;
    m21 *= invSx;
    m22 *= invSy;
    m23 *= invSz;
    m31 *= invSx;
    m32 *= invSy;
    m33 *= invSz;

    trace = m11 + m22 + m33;

    if (trace > 0) {
        s = 0.5 / mathf.sqrt(trace + 1);

        w = 0.25 / s;
        x = (m32 - m23) * s;
        y = (m13 - m31) * s;
        z = (m21 - m12) * s;
    } else if (m11 > m22 && m11 > m33) {
        s = 2 * mathf.sqrt(1 + m11 - m22 - m33);

        w = (m32 - m23) / s;
        x = 0.25 * s;
        y = (m12 + m21) / s;
        z = (m13 + m31) / s;
    } else if (m22 > m33) {
        s = 2 * mathf.sqrt(1 + m22 - m11 - m33);

        w = (m13 - m31) / s;
        x = (m12 + m21) / s;
        y = 0.25 * s;
        z = (m23 + m32) / s;
    } else {
        s = 2 * mathf.sqrt(1 + m33 - m11 - m22);

        w = (m21 - m12) / s;
        x = (m13 + m31) / s;
        y = (m23 + m32) / s;
        z = 0.25 * s;
    }

    rotation[0] = x;
    rotation[1] = y;
    rotation[2] = w;
    rotation[3] = z;

    return out;
};

mat4.setPosition = function(out, v) {
    var z = v[2];

    out[12] = v[0];
    out[13] = v[1];
    out[14] = z !== undefined ? z : 0;

    return out;
};

mat4.extractPosition = function(out, a) {

    out[12] = a[12];
    out[13] = a[13];
    out[14] = a[14];

    return out;
};

mat4.extractRotation = function(out, a) {
    var lx = vec3.lengthSqValues(a[0], a[1], a[2]),
        ly = vec3.lengthSqValues(a[4], a[5], a[6]),
        lz = vec3.lengthSqValues(a[8], a[9], a[10]),

        scaleX = lx !== 0 ? 1 / mathf.sqrt(lx) : lx,
        scaleY = ly !== 0 ? 1 / mathf.sqrt(ly) : ly,
        scaleZ = lz !== 0 ? 1 / mathf.sqrt(lz) : lz;

    out[0] = me[0] * scaleX;
    out[1] = me[1] * scaleX;
    out[2] = me[2] * scaleX;

    out[4] = me[4] * scaleY;
    out[5] = me[5] * scaleY;
    out[6] = me[6] * scaleY;

    out[8] = me[8] * scaleZ;
    out[9] = me[9] * scaleZ;
    out[10] = me[10] * scaleZ;

    return out;
};

mat4.extractRotationScale = function(out, a) {

    out[0] = a[0];
    out[1] = a[1];
    out[2] = a[2];

    out[4] = a[4];
    out[5] = a[5];
    out[6] = a[6];

    out[8] = a[8];
    out[9] = a[9];
    out[10] = a[10];

    return out;
};

mat4.translate = function(out, a, v) {
    var x = v[0],
        y = v[1],
        z = v[2],
        a00, a01, a02, a03,
        a10, a11, a12, a13,
        a20, a21, a22, a23;

    if (a === out) {
        out[12] = a[0] * x + a[4] * y + a[8] * z + a[12];
        out[13] = a[1] * x + a[5] * y + a[9] * z + a[13];
        out[14] = a[2] * x + a[6] * y + a[10] * z + a[14];
        out[15] = a[3] * x + a[7] * y + a[11] * z + a[15];
    } else {
        a00 = a[0];
        a01 = a[1];
        a02 = a[2];
        a03 = a[3];
        a10 = a[4];
        a11 = a[5];
        a12 = a[6];
        a13 = a[7];
        a20 = a[8];
        a21 = a[9];
        a22 = a[10];
        a23 = a[11];

        out[0] = a00;
        out[1] = a01;
        out[2] = a02;
        out[3] = a03;
        out[4] = a10;
        out[5] = a11;
        out[6] = a12;
        out[7] = a13;
        out[8] = a20;
        out[9] = a21;
        out[10] = a22;
        out[11] = a23;

        out[12] = a00 * x + a10 * y + a20 * z + a[12];
        out[13] = a01 * x + a11 * y + a21 * z + a[13];
        out[14] = a02 * x + a12 * y + a22 * z + a[14];
        out[15] = a03 * x + a13 * y + a23 * z + a[15];
    }

    return out;
};

mat4.scale = function(out, a, v) {
    var x = v[0],
        y = v[1],
        z = v[2];

    out[0] = a[0] * x;
    out[1] = a[1] * x;
    out[2] = a[2] * x;
    out[3] = a[3] * x;
    out[4] = a[4] * y;
    out[5] = a[5] * y;
    out[6] = a[6] * y;
    out[7] = a[7] * y;
    out[8] = a[8] * z;
    out[9] = a[9] * z;
    out[10] = a[10] * z;
    out[11] = a[11] * z;
    out[12] = a[12];
    out[13] = a[13];
    out[14] = a[14];
    out[15] = a[15];

    return out;
};

mat4.rotateX = function(out, a, angle) {
    var m12 = a[4],
        m22 = a[5],
        m32 = a[6],
        m42 = a[7],
        m13 = a[8],
        m23 = a[9],
        m33 = a[10],
        m43 = a[11],
        c = mathf.cos(angle),
        s = mathf.sin(angle);

    if (a !== out) {
        out[0] = a[0];
        out[1] = a[1];
        out[2] = a[2];
        out[3] = a[3];
        out[12] = a[12];
        out[13] = a[13];
        out[14] = a[14];
        out[15] = a[15];
    }

    out[4] = c * m12 + s * m13;
    out[5] = c * m22 + s * m23;
    out[6] = c * m32 + s * m33;
    out[7] = c * m42 + s * m43;

    out[8] = c * m13 - s * m12;
    out[9] = c * m23 - s * m22;
    out[10] = c * m33 - s * m32;
    out[11] = c * m43 - s * m42;

    return this;
};

mat4.rotateY = function(out, a, angle) {
    var m11 = a[0],
        m21 = a[1],
        m31 = a[2],
        m41 = a[3],
        m13 = a[8],
        m23 = a[9],
        m33 = a[10],
        m43 = a[11],
        c = mathf.cos(angle),
        s = mathf.sin(angle);

    if (a !== out) {
        out[4] = a[4];
        out[5] = a[5];
        out[6] = a[6];
        out[7] = a[7];
        out[12] = a[12];
        out[13] = a[13];
        out[14] = a[14];
        out[15] = a[15];
    }

    out[0] = c * m11 - s * m13;
    out[1] = c * m21 - s * m23;
    out[2] = c * m31 - s * m33;
    out[3] = c * m41 - s * m43;

    out[8] = c * m13 + s * m11;
    out[9] = c * m23 + s * m21;
    out[10] = c * m33 + s * m31;
    out[11] = c * m43 + s * m41;

    return this;
};

mat4.rotateZ = function(out, a, angle) {
    var m11 = a[0],
        m21 = a[1],
        m31 = a[2],
        m41 = a[3],
        m12 = a[4],
        m22 = a[5],
        m32 = a[6],
        m42 = a[7],
        c = mathf.cos(angle),
        s = mathf.sin(angle);

    if (a !== out) {
        out[8] = a[8];
        out[9] = a[9];
        out[10] = a[10];
        out[11] = a[11];
        out[12] = a[12];
        out[13] = a[13];
        out[14] = a[14];
        out[15] = a[15];
    }

    out[0] = c * m11 + s * m12;
    out[1] = c * m21 + s * m22;
    out[2] = c * m31 + s * m32;
    out[3] = c * m41 + s * m42;

    out[4] = c * m12 - s * m11;
    out[5] = c * m22 - s * m21;
    out[6] = c * m32 - s * m31;
    out[7] = c * m42 - s * m41;

    return out;
};

mat4.makeTranslation = function(out, v) {

    return mat4.set(
        out,
        1, 0, 0, v[0],
        0, 1, 0, v[1],
        0, 0, 1, v[2],
        0, 0, 0, 1
    );
};

mat4.makeScale = function(out, v) {

    return mat4.set(
        out,
        v[0], 0, 0, 0,
        0, v[1], 0, 0,
        0, 0, v[2], 0,
        0, 0, 0, 1
    );
};

mat4.makeRotationX = function(out, angle) {
    var c = mathf.cos(angle),
        s = mathf.sin(angle);

    return mat4.set(
        out,
        1, 0, 0, 0,
        0, c, -s, 0,
        0, s, c, 0,
        0, 0, 0, 1
    );
};

mat4.makeRotationY = function(out, angle) {
    var c = mathf.cos(angle),
        s = mathf.sin(angle);

    return mat4.set(
        out,
        c, 0, s, 0,
        0, 1, 0, 0, -s, 0, c, 0,
        0, 0, 0, 1
    );
};

mat4.makeRotationZ = function(out, angle) {
    var c = mathf.cos(angle),
        s = mathf.sin(angle);

    return mat4.set(
        out,
        c, -s, 0, 0,
        s, c, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    );
};

mat4.fromQuat = function(out, q) {
    var x = q[0],
        y = q[1],
        z = q[2],
        w = q[3],
        x2 = x + x,
        y2 = y + y,
        z2 = z + z,
        xx = x * x2,
        xy = x * y2,
        xz = x * z2,
        yy = y * y2,
        yz = y * z2,
        zz = z * z2,
        wx = w * x2,
        wy = w * y2,
        wz = w * z2;

    out[0] = 1 - (yy + zz);
    out[4] = xy - wz;
    out[8] = xz + wy;

    out[1] = xy + wz;
    out[5] = 1 - (xx + zz);
    out[9] = yz - wx;

    out[2] = xz - wy;
    out[6] = yz + wx;
    out[10] = 1 - (xx + yy);

    out[3] = 0;
    out[7] = 0;
    out[11] = 0;

    out[12] = 0;
    out[13] = 0;
    out[14] = 0;
    out[15] = 1;

    return out;
};

mat4.frustum = function(out, left, right, top, bottom, near, far) {
    var x = 2 * near / (right - left),
        y = 2 * near / (top - bottom),

        a = (right + left) / (right - left),
        b = (top + bottom) / (top - bottom),
        c = -(far + near) / (far - near),
        d = -2 * far * near / (far - near);

    out[0] = x;
    out[4] = 0;
    out[8] = a;
    out[12] = 0;
    out[1] = 0;
    out[5] = y;
    out[9] = b;
    out[13] = 0;
    out[2] = 0;
    out[6] = 0;
    out[10] = c;
    out[14] = d;
    out[3] = 0;
    out[7] = 0;
    out[11] = -1;
    out[15] = 0;

    return out;
};

mat4.perspective = function(out, fov, aspect, near, far) {
    var ymax = near * mathf.tan(fov * 0.5),
        ymin = -ymax,
        xmin = ymin * aspect,
        xmax = ymax * aspect;

    return mat4.frustum(out, xmin, xmax, ymax, ymin, near, far);
};

mat4.orthographic = function(out, left, right, top, bottom, near, far) {
    var w = right - left,
        h = top - bottom,
        p = far - near,

        x = (right + left) / w,
        y = (top + bottom) / h,
        z = (far + near) / p;

    out[0] = 2 / w;
    out[1] = 0;
    out[2] = 0;
    out[3] = 0;
    out[4] = 0;
    out[5] = 2 / h;
    out[6] = 0;
    out[7] = 0;
    out[8] = 0;
    out[9] = 0;
    out[10] = -2 / p;
    out[11] = 0;
    out[12] = -x;
    out[13] = -y;
    out[14] = -z;
    out[15] = 1;

    return out;
};

mat4.equal = function(a, b) {
    return !(
        a[0] !== b[0] ||
        a[1] !== b[1] ||
        a[2] !== b[2] ||
        a[3] !== b[3] ||
        a[4] !== b[4] ||
        a[5] !== b[5] ||
        a[6] !== b[6] ||
        a[7] !== b[7] ||
        a[8] !== b[8] ||
        a[9] !== b[9] ||
        a[10] !== b[10] ||
        a[11] !== b[11] ||
        a[12] !== b[12] ||
        a[13] !== b[13] ||
        a[14] !== b[14] ||
        a[15] !== b[15]
    );
};

mat4.notEqual = function(a, b) {
    return (
        a[0] !== b[0] ||
        a[1] !== b[1] ||
        a[2] !== b[2] ||
        a[3] !== b[3] ||
        a[4] !== b[4] ||
        a[5] !== b[5] ||
        a[6] !== b[6] ||
        a[7] !== b[7] ||
        a[8] !== b[8] ||
        a[9] !== b[9] ||
        a[10] !== b[10] ||
        a[11] !== b[11] ||
        a[12] !== b[12] ||
        a[13] !== b[13] ||
        a[14] !== b[14] ||
        a[15] !== b[15]
    );
};

mat4.str = function(out) {
    return (
        "Mat4[" + out[0] + ", " + out[4] + ", " + out[8] + ", " + out[12] + "]\n" +
        "     [" + out[1] + ", " + out[5] + ", " + out[9] + ", " + out[13] + "]\n" +
        "     [" + out[2] + ", " + out[6] + ", " + out[10] + ", " + out[14] + "]\n" +
        "     [" + out[3] + ", " + out[7] + ", " + out[11] + ", " + out[15] + "]"
    );
};
