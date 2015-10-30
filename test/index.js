var tape = require("tape"),
    mat4 = require("..");


tape("mat4.equal(a, b)", function(assert) {
    assert.equals(mat4.equal(mat4.create(), mat4.create()), true);
    assert.end();
});
