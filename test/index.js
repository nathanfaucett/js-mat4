var tape = require("tape"),
    mat4 = require("..");


tape("mat4.equal(a, b)", function(assert) {
    assert.equals(mat4.equals(mat4.create(), mat4.create()), true);
    assert.end();
});

tape("mat4.setRotationZ(out, a, angle)", function(assert) {
    assert.equals(mat4.equals(
        mat4.setRotationZ(mat4.create(), Math.PI / 2), mat4.create(
            0, -1, 0, 0,
            1, 0, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
        )
    ), true);
    assert.end();
});
