
// function that takes two arrays and createas an array of objects to use in chartjs plots
function convertData2ChartjsData(xin,yin) {
  var data = [];
  // var temp = {x: 987 , y: yin[0]};
  // console.log(xin)
  // console.log(yin)
  // console.log(xin.length)
  for (i = 0; i < xin.length; i++) {


    //    temp.x = xin[i];
    //    temp.y=yin[i];
    data.push({x:xin[i] , y:yin[i]});

    //     console.log(" i is equal to " + i)
    //     console.log(" temp is equal to " + temp.x)
    //     for (j = 0; j < i; j++) {
    //       console.log(" data is equal to " + data[j].x)
    //     };
    //     console.log(" data is equal to " + data[i].x)
    //
  };

  return data;                // Function returns data
};
