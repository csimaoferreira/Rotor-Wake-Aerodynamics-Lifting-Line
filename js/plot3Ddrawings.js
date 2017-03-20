


var camera, controls, scene, renderer;

// once everything is loaded, we run our Three.js stuff.
function maketestplot(width,height) {
    // create a scene, that will hold all our elements such as objects, cameras and lights.
    scene = new THREE.Scene();
    // create a camera, which defines where we're looking at.
    camera = new THREE.PerspectiveCamera(25, width / height, 0.3, 100000);
    // console.log(window);
    // create a render and set the size
    renderer = new THREE.WebGLRenderer();
    renderer.setClearColor(new THREE.Color(0xEEEEEE));
    renderer.setSize(width, height);
    var axes = new THREE.AxisHelper( 20 );
    scene.add(axes);

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
    // scene.add(sphere);


    var line= new createWakeMesh();
    scene.add(line);

    // position and point the camera to the center of the scene
    var csl=0.09;
    camera.position.x = -10/csl;
    camera.position.y = 10/csl;
    camera.position.z = 30/csl;

    camera.position.x = -0/csl;
    camera.position.y = 1/csl;
    camera.position.z = 30/csl;


    camera.lookAt(scene.position);


    // add the output of the renderer to the html element
    $("#WebGL-output").append(renderer.domElement);

    controls = new THREE.OrbitControls( camera, renderer.domElement );
    controls.addEventListener( 'change', render );
		// controls.enableZoom = false;

    // render the scene
    renderer.render(scene, camera);
};

function render() {

  renderer.render( scene, camera );

}


function createWakeMesh() {
  var s_Array = createArraySequence(0.,Math.PI/30, Math.PI);
  for (var i = 0; i < s_Array.length; i++) {
    s_Array[i]= (-1*(math.cos(s_Array[i])-1)/2*0.8+0.2)*(50);
  }
  // s_Array = createArraySequence(0.,1,100);


  var maxradius = math.max(s_Array);

  var theta_Array = createArraySequence(0.,Math.PI/50, 50*Math.PI);

  var rotor_wake_system = create_rotor_geometry(s_Array, maxradius, 6, 1, theta_Array,3);
  rings = rotor_wake_system.rings;

  var lineall = new THREE.Object3D();
 
  var geometry = new THREE.Geometry();
  var material = new THREE.LineBasicMaterial({
	color: 0xff0000
});

  for (var i = 0; i < rings.length; i++) {
   var filaments = rings[i].filaments;
   var nelements =(filaments.length-1)/2;

   var material = new THREE.LineBasicMaterial({
   	color: 0x0000ff
   });

  //  bound vortices
   var geometry = new THREE.Geometry();
   geometry.vertices.push(
   	new THREE.Vector3( filaments[0].x1, filaments[0].y1, filaments[0].z1 ),
    new THREE.Vector3( filaments[0].x2, filaments[0].y2, filaments[0].z2 )
   );
   var line = new THREE.Line( geometry, material );
   lineall.add(line);

 // trailing vortices from x1 of bound vortex
   var geometry1 = new THREE.Geometry();
      for (var  j= 1; j < 1+nelements; j++) {
    geometry1.vertices.push(
    	new THREE.Vector3( filaments[j].x2, filaments[j].y2, filaments[j].z2 )
    );
   }
   geometry1.vertices.push(
     new THREE.Vector3( filaments[nelements].x1, filaments[nelements].y1, filaments[nelements].z1 )
   );
       var line = new THREE.Line( geometry1, material );
    lineall.add(line)

    // trailing vortices from x2 of bound vortex
      var geometry1 = new THREE.Geometry();
         for (var  j= 1+nelements; j < filaments.length; j++) {
       geometry1.vertices.push(
       	new THREE.Vector3( filaments[j].x1, filaments[j].y1, filaments[j].z1 )
       );
      }
      geometry1.vertices.push(
        new THREE.Vector3( filaments[filaments.length-1].x2, filaments[filaments.length-1].y2, filaments[filaments.length-1].z2 )
      );
          var line = new THREE.Line( geometry1, material );
       lineall.add(line)


  };

// console.log(line);
  return(lineall)
}
