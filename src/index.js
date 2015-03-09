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
    var a00 = a[0],
        a01 = a[1],
        a02 = a[2],
        a03 = a[3],
        a10 = a[4],
        a11 = a[5],
        a12 = a[6],
        a13 = a[7],
        a20 = a[8],
        a21 = a[9],
        a22 = a[10],
        a23 = a[11],
        a30 = a[12],
        a31 = a[13],
        a32 = a[14],
        a33 = a[15],
        b0, b1, b2, b3;

    b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
    out[0] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
    out[1] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
    out[2] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
    out[3] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;

    b0 = b[4];
    b1 = b[5];
    b2 = b[6];
    b3 = b[7];
    out[4] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
    out[5] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
    out[6] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
    out[7] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;

    b0 = b[8];
    b1 = b[9];
    b2 = b[10];
    b3 = b[11];
    out[8] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
    out[9] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
    out[10] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
    out[11] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;

    b0 = b[12];
    b1 = b[13];
    b2 = b[14];
    b3 = b[15];
    out[12] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
    out[13] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
    out[14] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
    out[15] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;

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

mat4.determinant = function(a) {
    var a00 = a[0],
        a01 = a[1],
        a02 = a[2],
        a03 = a[3],
        a10 = a[4],
        a11 = a[5],
        a12 = a[6],
        a13 = a[7],
        a20 = a[8],
        a21 = a[9],
        a22 = a[10],
        a23 = a[11],
        a30 = a[12],
        a31 = a[13],
        a32 = a[14],
        a33 = a[15],

        b00 = a00 * a11 - a01 * a10,
        b01 = a00 * a12 - a02 * a10,
        b02 = a00 * a13 - a03 * a10,
        b03 = a01 * a12 - a02 * a11,
        b04 = a01 * a13 - a03 * a11,
        b05 = a02 * a13 - a03 * a12,
        b06 = a20 * a31 - a21 * a30,
        b07 = a20 * a32 - a22 * a30,
        b08 = a20 * a33 - a23 * a30,
        b09 = a21 * a32 - a22 * a31,
        b10 = a21 * a33 - a23 * a31,
        b11 = a22 * a33 - a23 * a32;

    return b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
};

mat4.inverse = function(out, a) {
    var a00 = a[0],
        a01 = a[1],
        a02 = a[2],
        a03 = a[3],
        a10 = a[4],
        a11 = a[5],
        a12 = a[6],
        a13 = a[7],
        a20 = a[8],
        a21 = a[9],
        a22 = a[10],
        a23 = a[11],
        a30 = a[12],
        a31 = a[13],
        a32 = a[14],
        a33 = a[15],

        b00 = a00 * a11 - a01 * a10,
        b01 = a00 * a12 - a02 * a10,
        b02 = a00 * a13 - a03 * a10,
        b03 = a01 * a12 - a02 * a11,
        b04 = a01 * a13 - a03 * a11,
        b05 = a02 * a13 - a03 * a12,
        b06 = a20 * a31 - a21 * a30,
        b07 = a20 * a32 - a22 * a30,
        b08 = a20 * a33 - a23 * a30,
        b09 = a21 * a32 - a22 * a31,
        b10 = a21 * a33 - a23 * a31,
        b11 = a22 * a33 - a23 * a32,

        det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;

    if (!det) {
        return mat4.identity(out);
    }
    det = 1.0 / det;

    out[0] = (a11 * b11 - a12 * b10 + a13 * b09) * det;
    out[1] = (a02 * b10 - a01 * b11 - a03 * b09) * det;
    out[2] = (a31 * b05 - a32 * b04 + a33 * b03) * det;
    out[3] = (a22 * b04 - a21 * b05 - a23 * b03) * det;
    out[4] = (a12 * b08 - a10 * b11 - a13 * b07) * det;
    out[5] = (a00 * b11 - a02 * b08 + a03 * b07) * det;
    out[6] = (a32 * b02 - a30 * b05 - a33 * b01) * det;
    out[7] = (a20 * b05 - a22 * b02 + a23 * b01) * det;
    out[8] = (a10 * b10 - a11 * b08 + a13 * b06) * det;
    out[9] = (a01 * b08 - a00 * b10 - a03 * b06) * det;
    out[10] = (a30 * b04 - a31 * b02 + a33 * b00) * det;
    out[11] = (a21 * b02 - a20 * b04 - a23 * b00) * det;
    out[12] = (a11 * b07 - a10 * b09 - a12 * b06) * det;
    out[13] = (a00 * b09 - a01 * b07 + a02 * b06) * det;
    out[14] = (a31 * b01 - a30 * b03 - a32 * b00) * det;
    out[15] = (a20 * b03 - a21 * b01 + a22 * b00) * det;

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

    out[0] = a[0] * scaleX;
    out[1] = a[1] * scaleX;
    out[2] = a[2] * scaleX;

    out[4] = a[4] * scaleY;
    out[5] = a[5] * scaleY;
    out[6] = a[6] * scaleY;

    out[8] = a[8] * scaleZ;
    out[9] = a[9] * scaleZ;
    out[10] = a[10] * scaleZ;

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
    var s = Math.sin(angle),
        c = Math.cos(angle),
        a10 = a[4],
        a11 = a[5],
        a12 = a[6],
        a13 = a[7],
        a20 = a[8],
        a21 = a[9],
        a22 = a[10],
        a23 = a[11];

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

    out[4] = a10 * c + a20 * s;
    out[5] = a11 * c + a21 * s;
    out[6] = a12 * c + a22 * s;
    out[7] = a13 * c + a23 * s;
    out[8] = a20 * c - a10 * s;
    out[9] = a21 * c - a11 * s;
    out[10] = a22 * c - a12 * s;
    out[11] = a23 * c - a13 * s;

    return out;
};

mat4.rotateY = function(out, a, angle) {
    var s = mathf.sin(angle),
        c = mathf.cos(angle),
        a00 = a[0],
        a01 = a[1],
        a02 = a[2],
        a03 = a[3],
        a20 = a[8],
        a21 = a[9],
        a22 = a[10],
        a23 = a[11];

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

    out[0] = a00 * c - a20 * s;
    out[1] = a01 * c - a21 * s;
    out[2] = a02 * c - a22 * s;
    out[3] = a03 * c - a23 * s;
    out[8] = a00 * s + a20 * c;
    out[9] = a01 * s + a21 * c;
    out[10] = a02 * s + a22 * c;
    out[11] = a03 * s + a23 * c;

    return out;
};

mat4.rotateZ = function(out, a, angle) {
    var s = mathf.sin(angle),
        c = mathf.cos(angle),
        a00 = a[0],
        a01 = a[1],
        a02 = a[2],
        a03 = a[3],
        a10 = a[4],
        a11 = a[5],
        a12 = a[6],
        a13 = a[7];

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

    out[0] = a00 * c + a10 * s;
    out[1] = a01 * c + a11 * s;
    out[2] = a02 * c + a12 * s;
    out[3] = a03 * c + a13 * s;
    out[4] = a10 * c - a00 * s;
    out[5] = a11 * c - a01 * s;
    out[6] = a12 * c - a02 * s;
    out[7] = a13 * c - a03 * s;

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
