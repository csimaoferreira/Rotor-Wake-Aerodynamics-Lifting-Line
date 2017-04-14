function solvepBEMLiftingLine(TSR,NELEMENTS, Nrotations) {
  // Nrotations=0.12
  // document.getElementById("message_calculate_circulation_BEM1").innerHTML = "Tip speed ratio = " +TSR.toFixed(2);
  var s_Array = createArraySequence(0.,(Math.PI)/NELEMENTS, Math.PI);
  // console.log(s_Array);
  // console.log(s_Array[1]+s_Array[s_Array.length-1]);
  var r_Array=[];
  for (var i = 0; i < s_Array.length; i++) {
    r_Array.push(-1*(math.cos(s_Array[i])-1)/2*0.8+0.2); // discretization for BEM model
  }

  for (var i = 0; i < s_Array.length; i++) {
    s_Array[i]= (-1*(math.cos(s_Array[i])-1)/2*0.8+0.2)*(50);
  }
  var maxradius = math.max(s_Array);
  // console.log(Nrotations);
  var theta_Array = createArraySequence(0.,Math.PI/10, Nrotations*2*Math.PI);
  var rotor_wake_system = create_rotor_geometry(s_Array, maxradius, TSR/(1-0.2), 1, theta_Array,3);
  // maketestplot(900, 400, rotor_wake_system,'WebGL-output')
  // console.log("LIFTING LINE");
  // console.log("LIFTING LINE");
  // console.log("LIFTING LINE");
  // console.log("LIFTING LINE");
  var resultsLL = solve_lifting_line_system_matrix_approach(rotor_wake_system, [1, 0, 0], TSR/50, maxradius);
  // console.log("BEM");
  // console.log("BEM");
  // console.log("BEM");
  // console.log("BEM");
  // console.log("BEM");
  // console.log("BEM");
  // console.log("BEM");
  // console.log("BEM");
  // console.log("BEM");
  var results = solveBEMmodel(1, r_Array, TSR*0.02 , 50, 3);
  var adim = Math.PI/(3*TSR/50);

  var results2 = calculateCT_CProtor_CPflow(results.a,results.aline,results.Fnorm, results.Ftan, 1, r_Array, TSR*0.02 , 50, 3);
  var CTCPliftline = calculateCT_CProtor_CPflow(resultsLL.a,resultsLL.aline,resultsLL.Fnorm, resultsLL.Ftan, 1, r_Array, TSR*0.02 , 50, 3);
  return({BEM: results , BEMloads: results2, LiftLine: resultsLL, LiftLineloads: CTCPliftline, rotor_wake_system: rotor_wake_system})

  // plotGamma2(results.r_R, nondim(results.Gamma, adim),  'BEM' , results.r_R, nondim(resultsLL.Gamma, adim),  'LiftLine' , "chart_calculate_circulation_Lift_BEM1")
  // document.getElementById("message_calc_circulation_BEM1_CTrotor").innerHTML =  results2.CTrotor.toFixed(2);
  // document.getElementById("message_calc_circulation_BEM1_CProtor").innerHTML =  results2.CProtor.toFixed(2);
  // document.getElementById("message_calc_circulation_LiftLine1_CTrotor").innerHTML =  CTCPliftline.CTrotor.toFixed(2);
  // document.getElementById("message_calc_circulation_LiftLine1_CProtor").innerHTML =  CTCPliftline.CProtor.toFixed(2);
}




function solvewingLiftingLine(Nelements,AspectRatio,Alpha) {
  var s_Array = createArraySequence(0.,(Math.PI)/Nelements, Math.PI);

  for (var i = 0; i < s_Array.length; i++) {
    s_Array[i]= (-1*(math.cos(s_Array[i])-1)/2)*(AspectRatio);
  }
  var WakeLength= 1000*AspectRatio;
  var rotor_wake_system = create_straight_wing_geometry(s_Array, Alpha, WakeLength);

  var results = solve_wing_lifting_line_system_matrix_approach(rotor_wake_system, Alpha);

  return({Cl: results[0], span: results[1], rotor_wake_system: rotor_wake_system})

}









function create_straight_wing_geometry(span_array, Alpha, WakeLength){
  var temp1
  var temp2
  var filaments = [];
  var ring = [];
  var infinity =WakeLength;
  var controlpoints = [];


  for (var i = 0; i < span_array.length-1; i++) {
    controlpoints.push( {coordinates: [ 0 , (span_array[i]+span_array[i+1])/2   , 0 ] , chord: 1, normal: [ 0,  0, 1] , tangential: [1, 0, 0]   } );
    temp1= {x1: infinity , y1:span_array[i], z1: infinity*Math.sin(Alpha*Math.PI/180), x2: 1.25 , y2:span_array[i], z2:0, Gamma: 0  }   ;
    filaments.push(temp1);
    temp1= {x1: 1.25 , y1:span_array[i], z1:0, x2:0 , y2:span_array[i], z2:0, Gamma: 0  }   ;
    filaments.push(temp1);
    temp1= {x1:0 , y1:span_array[i], z1:0, x2:0 , y2:span_array[i+1], z2:0, Gamma: 0  }   ;
    filaments.push(temp1);
    temp1= {x1:0 , y1:span_array[i+1], z1:0, x2:1.25 , y2:span_array[i+1], z2:0, Gamma: 0  }   ;
    filaments.push(temp1);
    temp1= {x1:1.25 , y1:span_array[i+1], z1:0, x2:infinity , y2:span_array[i+1], z2:infinity*Math.sin(Alpha*Math.PI/180), Gamma: 0  }   ;
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
  var bladepanels = [];
  var geodef; var geodef2;
  var r; var angle;  var angle2 ; var dx; var dz; var dy; var dtheta; var xt; var yt; var zt;
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

      // panel of the blade section
      geodef = geoBlade(span_array[i]/radius);
      angle = geodef[1]*Math.PI/180;
      geodef2 = geoBlade(span_array[i+1]/radius);
      angle2 = geodef2[1]*Math.PI/180;

      temp1= {
        p1: [-0.25*geodef[0]*Math.sin(-1*angle) , span_array[i], 0.25*geodef[0]*Math.cos(angle)],
        p2: [-0.25*geodef2[0]*Math.sin(-1*angle2) , span_array[i+1], 0.25*geodef2[0]*Math.cos(angle2)],
        p3: [0.75*geodef2[0]*Math.sin(-1*angle2) , span_array[i+1], -0.75*geodef2[0]*Math.cos(angle2)],
        p4: [0.75*geodef[0]*Math.sin(-1*angle) , span_array[i], -0.75*geodef[0]*Math.cos(angle)]
      };
      temp1.p1 = [ 0 ,  temp1.p1[1]*cosrot -  temp1.p1[2]*sinrot , temp1.p1[1]*sinrot +  temp1.p1[2]*cosrot ];
      temp1.p2 = [ 0 ,  temp1.p2[1]*cosrot -  temp1.p2[2]*sinrot , temp1.p2[1]*sinrot +  temp1.p2[2]*cosrot ];
      temp1.p3 = [ 0 ,  temp1.p3[1]*cosrot -  temp1.p3[2]*sinrot , temp1.p3[1]*sinrot +  temp1.p3[2]*cosrot ];
      temp1.p4 = [ 0 ,  temp1.p4[1]*cosrot -  temp1.p4[2]*sinrot , temp1.p4[1]*sinrot +  temp1.p4[2]*cosrot ];

      bladepanels.push(temp1);
    };

  };
  return({controlpoints: controlpoints ,  rings: ring, bladepanels:bladepanels});
};




    function solve_lifting_line_system_matrix_approach(rotor_wake_system,wind, Omega, rotorradius) {
      // this codes solves a lifting line model of a horizontal axis rotor
      // as inputs, it takes
      //      rotor_wake_system: data strucure that contains the geometry of the horseshoe vortex rings,
      //                         and the control points at the blade
      //      wind: unperturbed wind velocity, also known as U_infinity
      //      Omega: rotational velocity of the rotor
      //      rotorradius: the radius of the rotor

      // get controlpoints data structure
      var controlpoints = rotor_wake_system.controlpoints;
      // get horseshoe vortex rings data structure
      var rings = rotor_wake_system.rings;
      //
      // initialize variables that we will use during the calculation
      var velocity_induced =[]; // velocity induced by a horse vortex ring at a control point
      var up = []; var vp = []; var wp = []; // components of the velocity induced by one horseshoe vortex ring
      var u = 0;  var v = 0;  var w = 0; // total velocity induced at one control point
      var radialposition; var azimdir; // radial position of the control point
      var alpha; // angle of attack
      var GammaNew=[]; // new estimate of bound circulation
      var Gamma=[]; // current solution of bound circulation
      for (var i = 0; i < controlpoints.length; i++) { GammaNew.push(0);}; // initialize as zeros
      var vel1; var vmag; var vaxial; var vazim; var temploads; // velocities and loads at controlpoint
      var MatrixU = new Array(); // matrix of induction, for velocity component in x-direction
      var MatrixV = new Array(); // matrix of induction, for velocity component in y-direction
      var MatrixW = new Array(); // matrix of induction, for velocity component in z-direction
      // output variables
      var a_temp = new Array(); // output vector for axial induction
      var aline_temp = new Array();  // output vector for azimuthal induction
      var r_R_temp = new Array();  // output vector for radial position
      var Fnorm_temp = new Array();  // output vector for axial force
      var Ftan_temp = new Array();  // output vector for tangential force
      var Gamma_temp = new Array();  // output vector for circulation

      // the variables below are to setup the maximum number of iterations and convergence criteria
      var Niterations =1200;
      var errorlimit = 0.01;
      var error = 1.0; var refererror;
      var ConvWeight =0.3;

      // initalize and calculate matrices for velocity induced by horseshoe vortex rings
      // two "for cicles", each line varying wind controlpoint "icp", each column varying with
      // horseshoe vortex ring "jring"
      for (var icp= 0; icp < controlpoints.length; icp++) {
        MatrixU[icp] = new Array(); // new line of matrix
        MatrixV[icp] = new Array(); // new line of matrix
        MatrixW[icp] = new Array(); // new line of matrix
        for (var jring = 0; jring < rings.length; jring++) {
          // set ring strenth to unity, to calculate velocity induced by horseshoe vortex ring "jring"
          // at controlpoint "icp"
          rings[jring] = update_Gamma_sinle_ring(rings[jring],1,1);
          velocity_induced = velocity_induced_single_ring(rings[jring], controlpoints[icp].coordinates);
          // add compnent of velocity per unit strenght of circulation to induction matrix
          MatrixU[icp][jring] = velocity_induced[0];
          MatrixV[icp][jring] = velocity_induced[1];
          MatrixW[icp][jring] = velocity_induced[2];
        };
      };

      // calculate solution through an iteratuve process
      for (var  kiter = 0; kiter < Niterations; kiter++) {

        for (var ig = 0; ig < GammaNew.length; ig++) {
          Gamma[ig] = GammaNew[ig]; //update current bound circulation with new estimate
        }

        // calculate velocity, circulation and loads at the controlpoints
        for (var icp= 0; icp < controlpoints.length; icp++) {
          // determine radial position of the controlpoint;
          radialposition = Math.sqrt(math.dot(controlpoints[icp].coordinates, controlpoints[icp].coordinates));
          u=0; v=0; w=0; // initialize velocity
          // multiply icp line of Matrix with vector of circulation Gamma to calculate velocity at controlpoint
          for (var jring = 0; jring < rings.length; jring++) {
            u = u + MatrixU[icp][jring]*Gamma[jring]; // axial component of velocity
            v= v + MatrixV[icp][jring]*Gamma[jring]; // y-component of velocity
            w= w + MatrixW[icp][jring]*Gamma[jring]; // z-component of velocity
          };
          // calculate total perceived velocity
          vrot = math.cross([-Omega, 0 , 0]  , controlpoints[icp].coordinates ); // rotational velocity
          vel1 = [wind[0]+ u + vrot[0], wind[1]+ v + vrot[1] , wind[2]+ w + vrot[2]]; // total perceived velocity at section
          // calculate azimuthal and axial velocity
          azimdir = math.cross([-1/radialposition, 0 , 0]  , controlpoints[icp].coordinates ); // rotational direction
          vazim = math.dot(azimdir , vel1); // azimuthal direction
          vaxial =  math.dot([1, 0, 0] , vel1); // axial velocity
          // calculate loads using blade element theory
          temploads = loadBladeElement(vaxial, vazim, radialposition/rotorradius);
          // new point of new estimate of circulation for the blade section
          GammaNew[icp] = temploads[2];
          // update output vector
          a_temp[icp] =(-(u + vrot[0])/wind[0]);
          aline_temp[icp] =(vazim/(radialposition*Omega)-1);
          r_R_temp[icp] =(radialposition/rotorradius);
          Fnorm_temp[icp] =(temploads[0]);
          Ftan_temp[icp] =(temploads[1]);
          Gamma_temp[icp] =(temploads[2]);
        }; // end loop control points

        // check convergence of solution
        refererror =math.max(math.abs(GammaNew));
        refererror =Math.max(refererror,0.001); // define scale of bound circulation
        error =math.max(math.abs(math.subtract(GammaNew, Gamma))); // difference betweeen iterations
        error= error/refererror; // relative error
        if (error < errorlimit) {
          // if error smaller than limit, stop iteration cycle
          kiter=Niterations;
        }

        // set new estimate of bound circulation
        for (var ig = 0; ig < GammaNew.length; ig++) {
          GammaNew[ig] = (1-ConvWeight)*Gamma[ig] + ConvWeight*GammaNew[ig];
        }
      }; // end iteration loop

      // output results of converged solution
      return({a: a_temp , aline: aline_temp, r_R: r_R_temp, Fnorm: Fnorm_temp, Ftan: Ftan_temp , Gamma: Gamma_temp});
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
  // function to calculate the velocity induced by a straight 3D vortex filament
  // with circulation GAMMA at a point VP1. The geometry of the vortex filament
  // is defined by its edges: the filaments start at XV1 and ends at XV2.
  // the input CORE defines a vortex core radius, inside which the velocity
  // is defined  as a solid body rotation.
  // The function is adapted from the algorithm presented in:
  //                Katz, Joseph, and Allen Plotkin. Low-speed aerodynamics.
  //                Vol. 13. Cambridge university press, 2001.

  // read coordinates that define the vortex filament
  var X1 =XV1[0]; var Y1 =XV1[1]; var Z1 =XV1[2]; // start point of vortex filament
  var X2 =XV2[0]; var Y2 =XV2[1]; var Z2 =XV2[2]; // end point of vortex filament
  // read coordinates of target point where the velocity is calculated
  var XP =XVP1[0]; var YP =XVP1[1]; var ZP =XVP1[2];
  // calculate geometric relations for integral of the velocity induced by filament
  var R1=Math.sqrt(Math.pow((XP-X1), 2) + Math.pow((YP-Y1), 2) + Math.pow((ZP-Z1), 2) );
  var R2=Math.sqrt( Math.pow((XP-X2), 2) + Math.pow((YP-Y2), 2) + Math.pow((ZP-Z2), 2) );
  var R1XR2_X=(YP-Y1)*(ZP-Z2)-(ZP-Z1)*(YP-Y2);
  var R1XR2_Y=-(XP-X1)*(ZP-Z2)+(ZP-Z1)*(XP-X2);
  var R1XR2_Z=(XP-X1)*(YP-Y2)-(YP-Y1)*(XP-X2);
  var R1XR_SQR=Math.pow(R1XR2_X, 2)+ Math.pow(R1XR2_Y, 2)+ Math.pow(R1XR2_Z, 2);
  var R0R1 = (X2-X1)*(XP-X1)+(Y2-Y1)*(YP-Y1)+(Z2-Z1)*(ZP-Z1);
  var R0R2 = (X2-X1)*(XP-X2)+(Y2-Y1)*(YP-Y2)+(Z2-Z1)*(ZP-Z2);
  // check if target point is in the vortex filament core,
  // and modify to solid body rotation
  if (R1XR_SQR<Math.pow(CORE,2)) {
    R1XR_SQR=Math.pow(CORE,2);
    // GAMMA = 0;
  };
  if (R1<CORE) {
    R1=CORE;
    // GAMMA = 0;
  };
  if (R2<CORE) {
    R2=CORE;
    // GAMMA = 0;
  };
  // determine scalar
  var K=GAMMA/4/Math.PI/R1XR_SQR*(R0R1/R1 -R0R2/R2 );
  // determine the three velocity components
  var U=K*R1XR2_X;
  var V=K*R1XR2_Y;
  var W=K*R1XR2_Z;
  // output results, vector with the three velocity components
  var results = [U, V, W];
  return(results) };







  function solve_wing_lifting_line_system_matrix_approach(rotor_wake_system, Alpha){
    var controlpoints = rotor_wake_system.controlpoints;
    var rings = rotor_wake_system.rings;
    // initialize variabless
    var velocity_induced =[];
    var up = []; var vp = []; var wp = [];
    var u = 0;  var v = 0;  var w = 0;
    var alpha1, cl1, vmag1;
    var GammaNew=[];
    var ClNew=[];
    var rNew=[];
    var Gamma=[]; for (var i = 0; i < controlpoints.length; i++) { GammaNew.push(0);};
    for (var i = 0; i < controlpoints.length; i++) { ClNew.push(0);};
    for (var i = 0; i < controlpoints.length; i++) { rNew.push(controlpoints[i].coordinates[1]);};
    var Niterations = 340;
    var MatrixU = [];
    var MatrixV = [];
    var MatrixW = [];
    var errorlimit = 0.01;
    var error = 1.0; var refererror;
    var ConvWeight =0.1;

    // initalize and calculate matrix
    for (var icp= 0; icp < controlpoints.length; icp++) {
      for (var jring = 0; jring < rings.length; jring++) {
        rings[jring] = update_Gamma_sinle_ring(rings[jring],1,1);
        velocity_induced = velocity_induced_single_ring(rings[jring], controlpoints[icp].coordinates);
        up.push(velocity_induced[0]*1);
        vp.push(velocity_induced[1]*1);
        wp.push(velocity_induced[2]*1);
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
        // calculate total perceived velocity

        vel1 = [1  + u ,  v , 1*Math.sin(Alpha*Math.PI/180)+ w ]; // total perceived velocity at section
        angle1= Math.atan(vel1[2]/vel1[0]);
        ClNew[icp]=2*Math.PI*Math.sin(angle1);
        vmag = Math.sqrt(math.dot(vel1 , vel1));
        GammaNew[icp] = 0.5*1*vmag*ClNew[icp];

      }; // end loop control points
      refererror =math.max(math.abs(GammaNew));
      refererror =Math.max(refererror,0.001);
      // var errorold = error;
      error =math.max(math.abs(math.subtract(GammaNew, Gamma)));
      // console.log("error absolute " + error);
      error= error/refererror;
      ConvWeight = Math.max((1-error)*0.3,0.1);

      if (error<errorlimit) {
        kiter=Niterations;
      }


      for (var ig = 0; ig < GammaNew.length; ig++) {
        GammaNew[ig] = (1-ConvWeight)*Gamma[ig] + ConvWeight*GammaNew[ig];
      }
    }; // end iteration loop

    return [ClNew, rNew];
  };

























  //
