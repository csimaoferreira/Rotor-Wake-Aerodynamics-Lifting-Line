

// once everything is loaded, we run our Three.js stuff.
function maketestplot(width,height) {
    // create a scene, that will hold all our elements such as objects, cameras and lights.
    var scene = new THREE.Scene();
    // create a camera, which defines where we're looking at.
    var camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 1000);
    console.log(window);
    // create a render and set the size
    var renderer = new THREE.WebGLRenderer();
    renderer.setClearColor(new THREE.Color(0xEEEEEE));
    renderer.setSize(width, height);
    var axes = new THREE.AxisHelper( 20 );
    // scene.add(axes);

    // create the ground plane
    var planeGeometry = new THREE.PlaneGeometry(60,20);
    var planeMaterial = new THREE.MeshBasicMaterial({color: 0xcccccc});
    var plane = new THREE.Mesh(planeGeometry,planeMaterial);
    // rotate and position the plane
    plane.rotation.x=-0.5*Math.PI;
    plane.position.x=15
    plane.position.y=0
    plane.position.z=0
    // add the plane to the scene
    // scene.add(plane);
    // create a cube
    var cubeGeometry = new THREE.CubeGeometry(4,4,4);
    var cubeMaterial = new THREE.MeshBasicMaterial({color: 0xff0000, wireframe: true});
    var cube = new THREE.Mesh(cubeGeometry, cubeMaterial);
    // position the cube
    cube.position.x=-4;
    cube.position.y=3;
    cube.position.z=0;
    // add the cube to the scene
    // scene.add(cube);
    var sphereGeometry = new THREE.SphereGeometry(4,20,20);
    var sphereMaterial = new THREE.MeshBasicMaterial({color: 0x7777ff, wireframe: true});
    var sphere = new THREE.Mesh(sphereGeometry,sphereMaterial);
    // position the sphere
    sphere.position.x=20;
    sphere.position.y=4;
    sphere.position.z=2;
    // add the sphere to the scene
    scene.add(sphere);


    var line= new createWakeMesh();
    scene.add(line);

    // position and point the camera to the center of the scene
    camera.position.x = -30;
    camera.position.y = 40;
    camera.position.z = 30;
    camera.lookAt(scene.position);
    // add the output of the renderer to the html element
    $("#WebGL-output").append(renderer.domElement);
    // render the scene
    renderer.render(scene, camera);
};



function createWakeMesh() {
  var s_Array = createArraySequence(0.,Math.PI/5, Math.PI);
  for (var i = 0; i < s_Array.length; i++) {
    s_Array[i]= -1*(math.cos(s_Array[i])-1)/2*(5);
  }
  // s_Array = createArraySequence(0.,1,100);

  var theta_Array = createArraySequence(0.,Math.PI/10, 5*Math.PI);

  var rotor_wake_system = create_rotor_geometry(s_Array, 5, 6, 1, theta_Array);
  rings = rotor_wake_system.rings;

  var lineall = new THREE.Object3D();

  var geometry = new THREE.Geometry();
  var material = new THREE.LineBasicMaterial({
	color: 0xff0000
});

  for (var i = 0; i < 1; i++) {
   var filaments = rings[i].filaments;
   for (var  j= 0; j < filaments.length; j++) {
     geometry.vertices.push(THREE.Vector3( filaments[j].x1, filaments[j].y1, filaments[j].z1) ,THREE.Vector3(filaments[j].x2, filaments[j].y2, filaments[j].z2) );
     var line = new THREE.Line( geometry, material );
     lineall.add(line)
     geometry = new THREE.Geometry();

   }
  }







  return(lineall)
}
