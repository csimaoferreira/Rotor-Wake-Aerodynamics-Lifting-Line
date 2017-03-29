


var camera, controls, scene, renderer;
var rotation_velocity = 2*3.14/5;
var pointerAnimationrequest;

// once everything is loaded, we run our Three.js stuff.
function maketestplot(width,height) {
    // create a scene, that will hold all our elements such as objects, cameras and lights.
    scene = new THREE.Scene();
    // create a camera, which defines where we're looking at.
    camera = new THREE.PerspectiveCamera(25, width / height, 0.3, 100000);
    // console.log(window);
    // create a render and set the size
    renderer = new THREE.WebGLRenderer();
    renderer.setClearColor(new THREE.Color(0xFFFFFF));
    renderer.setSize(width, height);

    var s_Array = createArraySequence(0.,Math.PI/10, Math.PI);
    for (var i = 0; i < s_Array.length; i++) {
      s_Array[i]= (-1*(math.cos(s_Array[i])-1)/2*0.8+0.2)*(50);
    }
    var maxradius = math.max(s_Array);
    var theta_Array = createArraySequence(0.,Math.PI/20, 50*Math.PI);
    var rotor_wake_system = create_rotor_geometry(s_Array, maxradius, 6, 1, theta_Array,3);
    rings = rotor_wake_system.rings;




    var rotor = createRotor(rotor_wake_system.bladepanels);

    var wake= new createWakeMesh(rings);
    // scene.add(wake);

    var rotorwake = new THREE.Object3D();
    rotorwake.name = 'RotorWake';
    rotorwake.add(wake);
    rotorwake.add(rotor);

    scene.add(rotorwake);

    // position and point the camera to the center of the scene
    var csl=0.09;
    // camera.position.x = -20/csl;
    // camera.position.y = 10/csl;
    // camera.position.z = 30/csl;

    camera.position.x = -100;
    camera.position.y = 10;
    camera.position.z = 300;
    scene.position.x= 50;

    camera.lookAt(scene.position);

    // add the output of the renderer to the html element
    $("#WebGL-output").append(renderer.domElement);

    // controls = new THREE.OrbitControls( camera, renderer.domElement );
    // controls.addEventListener( 'change', render );
		// controls.enableZoom = false;

    // render the scene
    renderer.render(scene, camera);
};

function render() {
  var d = new Date();
  var time = (d.getTime())/1000;
  var rotation = rotation_velocity*time;
  // var rotation = rotation_velocity*time;

  rotation = rotation % (2*Math.PI);
  var RotorWake = scene.getObjectByName( "RotorWake", true );
  RotorWake.rotation.x=rotation;

  // var idslide=Reveal.getCurrentSlide().id;
  // console.log(idslide);
  // console.log(idslide);
  // if (idslide=="Coverslide") {
  //   console.log("we are here");
  //   console.log("we are here");
  //   console.log(idslide);
  //   console.log(idslide);
    pointerAnimationrequest=requestAnimationFrame(render);
  // };


  renderer.render( scene, camera );
  // controls.update();

  // console.log("rotation is " + rotation);
}

// render();

function createWakeMesh(rings) {

  var lineall = new THREE.Object3D();
  lineall.name = 'WakeMesh'

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




function createRotor(bladesections) {

  var rotor = new THREE.Object3D();
  rotor.name = 'Rotor'

  var geometry = new THREE.Geometry();

  var material = new THREE.MeshBasicMaterial({
                       color:0xfed731,
                       side:THREE.DoubleSide
                   });

  for (var i = 0; i < bladesections.length; i++) {

    geometry.vertices.push(
    	new THREE.Vector3( bladesections[i].p1[0],  bladesections[i].p1[1], bladesections[i].p1[2] ),
      new THREE.Vector3( bladesections[i].p2[0],  bladesections[i].p2[1], bladesections[i].p2[2] ),
      new THREE.Vector3( bladesections[i].p3[0],  bladesections[i].p3[1], bladesections[i].p3[2] ),
      new THREE.Vector3( bladesections[i].p4[0],  bladesections[i].p4[1], bladesections[i].p4[2] )
    );

    geometry.faces.push( new THREE.Face3( i*4+0, i*4+1, i*4+2 ) );
    geometry.faces.push( new THREE.Face3( i*4+0, i*4+2, i*4+3 ) );
  };

  var Mesh = new THREE.Mesh(geometry, material);
  rotor.add(Mesh);

// console.log(line);
  return(rotor)
}
