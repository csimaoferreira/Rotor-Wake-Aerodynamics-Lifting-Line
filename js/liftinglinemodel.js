 
function create_straight_wing_geometry(span_array){
  var temp1
  var temp2
  var filaments = [];
  var ring = [];
  var infinity =10000000000000;
  var controlpoints = [];


  for (var i = 0; i < span_array.length-1; i++) {
    controlpoints.push( {coordinates: [ 0 , (span_array[i]+span_array[i+1])/2   , 0 ] , chord: 1, normal: [ 0,  0, 1] , tangential: [1, 0, 0]   } );
    temp1= {x1: infinity , y1:span_array[i], z1:0, x2:0 , y2:span_array[i], z2:0, Gamma: 0  }   ;
    filaments.push(temp1);
    temp1= {x1:0 , y1:span_array[i], z1:0, x2:0 , y2:span_array[i+1], z2:0, Gamma: 0  }   ;
    filaments.push(temp1);
    temp1= {x1:0 , y1:span_array[i+1], z1:0, x2:infinity , y2:span_array[i+1], z2:0, Gamma: 0  }   ;
    filaments.push(temp1);
    ring.push({filaments: filaments});
    filaments = [];
  };
  return({controlpoints: controlpoints ,  rings: ring});
};


function geoBlade(r_R) {
  var radial = [0, 0.3, .5, .8, 1];
  // var chorddist = [.05, .04, .03, .02, .015];
  // var twistdist = [-12, -8, -5, -4, 0.];
  var pitch = 2;
  var chord = 3*(1-r_R)+1;
  var twist = -14*(1-r_R);
  var result =[chord , twist + pitch];
  return result;
};


function create_rotor_geometry(span_array, radius, tipspeedratio, Uinf, theta_array, nblades) {
  // create rotor eometry and rotor-wake circulation system
  var filaments = [];
  var ring = [];
  var controlpoints = [];
  var geodef;
  var r; var angle; var dx; var dz; var dy; var dtheta; var xt; var yt; var zt;
  var temp1; var temp2;

  for (var krot = 0; krot < nblades; krot++) {

    var angle_rotation = 2*Math.PI/nblades*krot;
    var cosrot = Math.cos(angle_rotation);
    var sinrot = Math.sin(angle_rotation);

    for (var i = 0; i < span_array.length-1; i++) {
      r = (span_array[i]+span_array[i+1])/2;
      geodef = geoBlade(r/radius);
      angle = geodef[1]*Math.PI/180;
      // define controlpoints
      temp1 = {coordinates: [ 0 ,  r  , 0 ] , chord: geodef[0], normal: [ Math.cos(angle),  0, -1*Math.sin(angle)] , tangential: [-1*Math.sin(angle), 0, -1*Math.cos(angle)]   };
      // rotate blade to position
      temp1.coordinates = [ 0 ,  temp1.coordinates[1]*cosrot -  temp1.coordinates[2]*sinrot , temp1.coordinates[1]*sinrot +  temp1.coordinates[2]*cosrot ];
      temp1.normal = [ temp1.normal[0],  temp1.normal[1]*cosrot -  temp1.normal[2]*sinrot , temp1.normal[1]*sinrot +  temp1.normal[2]*cosrot ];
      temp1.tangential = [ temp1.tangential[0] ,  temp1.tangential[1]*cosrot -  temp1.tangential[2]*sinrot , temp1.tangential[1]*sinrot +  temp1.tangential[2]*cosrot ];

      controlpoints.push( temp1 );
      // define bound vortex filament
      temp1= {x1:0 , y1:span_array[i], z1:0, x2:0 , y2:span_array[i+1], z2:0, Gamma: 0  }   ;
      // rotate filament to position

      filaments.push(temp1);
      // create trailing filaments, at x1 of bound filament
      geodef = geoBlade(span_array[i]/radius);
      angle = geodef[1]*Math.PI/180;
      temp1= {x1: geodef[0]*Math.sin(-1*angle), y1:span_array[i], z1: -1*geodef[0]*Math.cos(angle), x2:0 , y2:span_array[i], z2:0, Gamma: 0  }   ;
      filaments.push(temp1);
      for (var j = 0; j < theta_array.length-1; j++) {
        xt = filaments[filaments.length-1].x1;
        yt = filaments[filaments.length-1].y1;
        zt = filaments[filaments.length-1].z1;
        dy = (Math.cos(-theta_array[j+1])-Math.cos(-theta_array[j])) * span_array[i];
        dz = (Math.sin(-theta_array[j+1])-Math.sin(-theta_array[j])) * span_array[i];
        dx = (theta_array[j+1]-theta_array[j])/tipspeedratio*radius;

        temp1= {x1: xt+dx, y1: yt+dy, z1: zt+dz, x2: xt , y2:yt, z2:zt, Gamma: 0  }   ;
        // rotate filament to position

        filaments.push(temp1);
      };

      // create trailing filaments, at x2 of bound filament
      geodef = geoBlade(span_array[i+1]/radius);
      angle = geodef[1]*Math.PI/180;
      temp1= {x1:0 , y1:span_array[i+1], z1: 0, x2:geodef[0]*Math.sin(-1*angle) , y2:span_array[i+1], z2:-1*geodef[0]*Math.cos(angle), Gamma: 0  }   ;
      filaments.push(temp1);
      for (var j = 0; j < theta_array.length-1; j++) {
        xt = filaments[filaments.length-1].x2;
        yt = filaments[filaments.length-1].y2;
        zt = filaments[filaments.length-1].z2;
        dy = (Math.cos(-theta_array[j+1])-Math.cos(-theta_array[j])) * span_array[i+1];
        dz = (Math.sin(-theta_array[j+1])-Math.sin(-theta_array[j])) * span_array[i+1];
        dx = (theta_array[j+1]-theta_array[j])/tipspeedratio*radius;

        temp1= {x1: xt, y1: yt, z1: zt, x2: xt+dx , y2:yt+dy, z2:zt+dz, Gamma: 0  }   ;
        // rotate filament to position

        filaments.push(temp1);
      };


      for (var ifil = 0; ifil < filaments.length; ifil++) {
        temp1=filaments[ifil];
        temp2 = [ temp1.y1*cosrot -  temp1.z1*sinrot , temp1.y1*sinrot +  temp1.z1*cosrot , temp1.y2*cosrot -  temp1.z2*sinrot , temp1.y2*sinrot +  temp1.z2*cosrot ];
        temp1.y1 = temp2[0];
        temp1.z1 = temp2[1];
        temp1.y2 = temp2[2];
        temp1.z2 = temp2[3];
        filaments[ifil] = temp1;
      }




      ring.push({filaments: filaments});
      filaments = [];
    };

  };
  return({controlpoints: controlpoints ,  rings: ring});
};


function solve_lifting_line_system_matrix_approach(rotor_wake_system,wind, Omega) {
  var controlpoints = rotor_wake_system.controlpoints;
  var rings = rotor_wake_system.rings;
  var velocity_induced =[];
  var up = [] ;
  var vp = [];
  var wp = [];
  var u = 0;
  var v = 0;
  var w = 0;
  var alpha;
  var GammaNew=[];
  var Gamma=[];
  for (var i = 0; i < controlpoints.length; i++) {
    GammaNew.push(0);
  };

  // console.log("gam1 is " + GammaNew);
  var vel1; var vmag; var vn ;
  var cl;
  var vtan; var vnorm; var vall; var vrot;
  var Niterations = 140;
  var MatrixU = [];
  var MatrixV = [];
  var MatrixW = [];
  var Niterations =1200;
  var errorlimit = 0.001;
  var error = 1.0; var refererror;
  var ConvWeight =0.1;

  // initalize and calculate matrix
  for (var icp= 0; icp < controlpoints.length; icp++) {
    for (var jring = 0; jring < rings.length; jring++) {
      rings[jring] = update_Gamma_sinle_ring(rings[jring],1,1);
      velocity_induced = velocity_induced_single_ring(rings[jring], controlpoints[icp].coordinates);
      // console.log(velocity_induced[0]);
      // var a = velocity_induced[0]*1
      // u = velocity_induced[0]*1;
      up.push(velocity_induced[0]*1);
      vp.push(velocity_induced[1]*1);
      wp.push(velocity_induced[2]*1);
      // console.log(wp);
      velocity_induced =[];
    };
      MatrixU.push(up);
      MatrixV.push(vp);
      MatrixW.push(wp);
      up =[]; vp =[]; wp =[];
  };

  for (var  kiter = 0; kiter < Niterations; kiter++) {
    Gamma=[];
    for (var ig = 0; ig < GammaNew.length; ig++) {
      Gamma.push(GammaNew[ig]);
    }
    for (var icp= 0; icp < controlpoints.length; icp++) {
      u=0; v=0; w=0;
      for (var jring = 0; jring < rings.length; jring++) {
        u = u + MatrixU[icp][jring]*Gamma[jring];
        v= v + MatrixV[icp][jring]*Gamma[jring];
        w= w + MatrixW[icp][jring]*Gamma[jring];
      };
      vrot = math.cross([-Omega, 0 , 0]  , controlpoints[icp].coordinates );
      vel1 = [wind[0]+ u + vrot[0], wind[1]+ v + vrot[1] , wind[2]+ w + vrot[2]];
      vtan = math.dot(controlpoints[icp].tangential , vel1);
      vnorm = math.dot(controlpoints[icp].normal , vel1);
      vmag = Math.sqrt(math.dot(vel1 , vel1));
      alpha = Math.atan(vnorm/vtan);
      cl =2*Math.PI*Math.sin(alpha);
      GammaNew[icp] = cl*0.5*vmag*controlpoints[icp].chord;
    }; // end loop control points
    refererror =math.max(math.abs(GammaNew));
    refererror =Math.max(refererror,0.001);
    // var errorold = error;
    error =math.max(math.abs(math.subtract(GammaNew, Gamma)));
    // console.log("error absolute " + error);
    error= error/refererror;

      var temp1 = (kiter+1)*1.0;
      ConvWeight = Math.max((1-error)/((Math.max(controlpoints.length/10 , ((kiter+1)/10.0))*2.0)) , 0.5/(Math.max(controlpoints.length , ((kiter+1)*1.0))));
      ConvWeight = Math.max((1-error)/(Math.sqrt(controlpoints.length*1)), 1/(Math.sqrt(controlpoints.length*1)));

      if (controlpoints.length<60) {
        ConvWeight = Math.max((1-error)*0.3,0.1);
      } else {
        ConvWeight =  Math.max((1-error),0.5)/Math.pow(controlpoints.length , 8/10);

      }
      // ConvWeight =  Math.max(0.5,0.5)/Math.pow(controlpoints.length , 7/10);
      // ConvWeight = 1/Math.pow(controlpoints.length , 1.1);
      // ConvWeight = Math.min(ConvWeight , 0.05);
      // ConvWeight = 0.5*(ConvWeight2 + ConvWeight);



    // console.log(error*100 + "%  in iteration " + kiter);
    if (error<errorlimit) {
      kiter=Niterations;
      console.log("very small error achieved");
    }

    console.log(error*100 + "%  in iteration " + kiter + " with conv " + ConvWeight);


    for (var ig = 0; ig < GammaNew.length; ig++) {
      GammaNew[ig] = (1-ConvWeight)*Gamma[ig] + ConvWeight*GammaNew[ig];
    }
// (1+Math.random()*0.4-0.2)*
  }; // end iteration loop
  console.log(GammaNew);


  return(GammaNew);
};





function update_Gamma_sinle_ring(ring,GammaNew,WeightNew) {

  for (var i = 0; i < ring.filaments.length; i++) {
    ring.filaments[i].Gamma = ring.filaments[i].Gamma*(1-WeightNew) + WeightNew * GammaNew;
  };
  return(ring);
};



function velocity_induced_rings(rings,controlpoints) {
  var tempvel1=[0 ,0, 0];
  var velind=[0 ,0, 0];
  for (var i = 0; i < rings.length; i++) {
    tempvel1 =velocity_induced_sinle_ring(rings[i],controlpoints);
    velind = math.add(velind,tempvel1);
  };
  // console.log(velind)

  return(velind);

};


function velocity_induced_single_ring(ring,controlpoint) {
  var tempvel1=[0 ,0, 0];
  var velind=[0 ,0, 0];
  var GAMMA;
  var XV1;
  var XV2;
  var XVP;
  var CORE = 0.00001;
  for (var i = 0; i < ring.filaments.length; i++) {
    GAMMA = ring.filaments[i].Gamma;
    XV1 = [ ring.filaments[i].x1 , ring.filaments[i].y1   ,  ring.filaments[i].z1 ];
    XV2 = [ ring.filaments[i].x2 , ring.filaments[i].y2   ,  ring.filaments[i].z2 ];
    XVP = controlpoint;
    tempvel1 =velocity_3D_from_vortex_filament(GAMMA,XV1, XV2, controlpoint ,CORE);
    //  console.log(i);
    //  console.log(tempvel1)
    //  console.log("XV1 " + XV1 + " XV2 " + XV2 + " XVP " + XVP );
    // //  console.log(velind)
    // velind = math.add(velind,tempvel1);
    velind[0]+=tempvel1[0];
    velind[1]+=tempvel1[1];
    velind[2]+=tempvel1[2];
      // console.log(velind)
  };
  return velind;
};


// 3D velocity induced by a vortex filament


function velocity_3D_from_vortex_filament(GAMMA,XV1, XV2, XVP1,CORE){
  var X1 =XV1[0];
  var Y1 =XV1[1];
  var Z1 =XV1[2];
  var X2 =XV2[0];
  var Y2 =XV2[1];
  var Z2 =XV2[2];
  var XP =XVP1[0];
  var YP =XVP1[1];
  var ZP =XVP1[2];

  var ERR=CORE;




  var R1=Math.sqrt(Math.pow((XP-X1), 2) + Math.pow((YP-Y1), 2) + Math.pow((ZP-Z1), 2) );
  var R2=Math.sqrt( Math.pow((XP-X2), 2) + Math.pow((YP-Y2), 2) + Math.pow((ZP-Z2), 2) );

  if (R1<ERR) {
    R1=ERR;
    GAMMA = 0;
  };

  if (R2<ERR) {
    R2=ERR;
    GAMMA = 0;

  };

  var R1XR2_X=(YP-Y1)*(ZP-Z2)-(ZP-Z1)*(YP-Y2);
  var R1XR2_Y=-(XP-X1)*(ZP-Z2)+(ZP-Z1)*(XP-X2);
  var R1XR2_Z=(XP-X1)*(YP-Y2)-(YP-Y1)*(XP-X2);

  var R1XR_SQR=Math.pow(R1XR2_X, 2)+ Math.pow(R1XR2_Y, 2)+ Math.pow(R1XR2_Z, 2);

  if (R1XR_SQR<Math.pow(ERR,2)) {
    R1XR_SQR=Math.pow(ERR,2);
    GAMMA = 0;

  };

  var R0R1 = (X2-X1)*(XP-X1)+(Y2-Y1)*(YP-Y1)+(Z2-Z1)*(ZP-Z1);
  var R0R2 = (X2-X1)*(XP-X2)+(Y2-Y1)*(YP-Y2)+(Z2-Z1)*(ZP-Z2);

  var K=GAMMA/4/Math.PI/R1XR_SQR*(R0R1/R1 -R0R2/R2 );

  var U=K*R1XR2_X;
  var V=K*R1XR2_Y;
  var W=K*R1XR2_Z;


// console.log("GAMMMA");
// console.log(GAMMA);
//
//   console.log("U V W");
//   console.log([U, V, W]);

  var results = [U, V, W];
  return(results) };





























  //
