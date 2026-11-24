// The mouse turns the patient in your hand instead of moving the camera: drag to turn it about its own middle,
// double-click to put it back. The wheel still zooms and the right button still pans. (The scene puts the middle of
// the body at the scene's origin, so the patient turns about itself.)
var T = S.THREE, pivot = new T.Group(), drag = null, el = S.renderer.domElement;
pivot.name = 'patient-turntable';
S.root.add(pivot);
S.controls.enableRotate = false;

function holder() {
  var h = S.root.getObjectByName('patient');
  if (h && h.parent !== pivot) pivot.add(h);     // (again after the scene is reloaded)
  return h;
}
function turn(axisWorld, angle) {
  var rq = S.root.getWorldQuaternion(new T.Quaternion()).invert();
  var axis = axisWorld.clone().applyQuaternion(rq).normalize();
  pivot.quaternion.premultiply(new T.Quaternion().setFromAxisAngle(axis, angle));
}
el.addEventListener('pointerdown', function (e) {
  if (e.button !== 0) return;
  drag = { x: e.clientX, y: e.clientY };
});
window.addEventListener('pointermove', function (e) {
  if (!drag) return;
  var dx = e.clientX - drag.x, dy = e.clientY - drag.y, k = Math.PI / Math.max(300, el.clientHeight);
  drag.x = e.clientX; drag.y = e.clientY;
  var cam = S.camera, up = new T.Vector3(0, 1, 0).applyQuaternion(cam.quaternion), right = new T.Vector3(1, 0, 0).applyQuaternion(cam.quaternion);
  turn(up, dx * k);
  turn(right, dy * k);
});
window.addEventListener('pointerup', function () { drag = null; });
el.addEventListener('dblclick', function () { pivot.quaternion.identity(); });
S.on('step', function () { holder(); });
