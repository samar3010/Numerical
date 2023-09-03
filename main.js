var resultText = document.getElementById("resultText");
var eps, maxIter;
var fx = "-0.6 * x^2 + 2.4 * x + 5.5";

function convert(){
  if(!fx) return;
  while(fx.includes("^")){
    fx = fx.replace("^", "**");
  }
}

// ----- CHAPTER 1 -------

function clearr(arr){
  arr.forEach(e => {
    document.getElementById(e).value = "";
  })
}

//Initializing desired error % & max iterations
function init(){
  if (document.getElementById("inputMethod").value != "") fx = document.getElementById("inputMethod").value;
  eps = parseFloat(document.getElementById("epsText").value) || 0;
  maxIter = parseInt(document.getElementById("iterText").value) || Infinity;
}
function f(x){
  let ret;
  if (typeof math == 'undefined') ret = roundTo(eval(fx), 3);
  else if(typeof math != 'undefined') ret = roundTo(math.eval(fx, {x: x}), 3);
  return ret;
}

function fDash(x){
  let xDash = fx;
  while(xDash.includes('**')) xDash = xDash.replace('**', '^');
  derv = math.derivative(xDash, 'x');
  return roundTo(derv.eval({x: x}), 3);
}

//Rounding (to third decimal point)
function r(x){
  if(isNaN(x)) return x;
  return roundTo(x, 3);
}

//Bisection Method
function bisect(xl, xu){
  init();
  let iterations = ""
  let xrOld = 0, error = 0, xr = 0, iter = 0;
  xl = parseFloat(xl);
  xu = parseFloat(xu);
  
  if(f(xl) * f(xu) > 0) return resultText.value = "Input error.";

  do{
      xrOld = xr;
      xr = r( (xl + xu)/2 );
      if(iter == 0) error = 100;
      else error = r( Math.abs((xr - xrOld)/xr) * 100 );
      iterations = `${iterations}-=Iteration: ${iter}=-\nXl: ${xl}\nf(xl): ${f(xl)}\nXu: ${xu}\nf(xu): ${f(xu)}\nxr: ${xr}\nf(xr): ${f(xr)}\nError%: ${error}%\n\n`;
      if (f(xl) * f(xr) > 0){
          xl = xr;
      }else xu = xr;
      iter++;
  }while (error > eps && iter < maxIter);
  iterations = `${iterations}\n\nRoot = ${xr}`;
  resultText.value = iterations;
}

//False Position Method
function falsePosition(xl, xu){
  init();
  let iterations = "";
  let xrOld = 0, error = 0, xr = 0, iter = 0;
  xl = parseFloat(xl);
  xu = parseFloat(xu);
  if(f(xl) * f(xu) > 0) return resultText.value = "Input error.";
  do{
      xrOld = xr;
      xr = r( xu - (f(xu) * (xl - xu)) / (f(xl) - f(xu)) );
      if(iter == 0) error = 100;
      else error = r( Math.abs((xr - xrOld)/xr) * 100 );
      iterations = `${iterations}-=Iteration: ${iter}=-\nXl: ${xl}\nf(xl): ${f(xl)}\nXu: ${xu}\nf(xu): ${f(xu)}\nxr: ${xr}\nf(xr): ${f(xr)}\nError%: ${error}%\n\n`;
      if (f(xl) * f(xr) > 0){
          xl = xr;
      }else xu = xr;
      iter++;
  }while (error > eps && iter < maxIter);
  iterations = `${iterations}\n\nRoot = ${xr}`;
  resultText.value = iterations;
}

//Simple Fixed Point Method
function simpleFixed(xi){
  init();
  let iterations = "";
  xi = parseFloat(xi), iter = 0, xiPlus1 = 0, error = 0;
  do{
    if(iter != 0) xi = xiPlus1;
    xiPlus1 = r( f(xi) );
    console.log(iterations);
    if(iter == 0) error = 100;
    else error = r( Math.abs((xiPlus1 - xi) / xiPlus1) * 100 );
    iterations = `${iterations}=Iteration: ${iter}=-\nXl: ${xi}\nf(xi): ${f(xi)}\nError%: ${error}%\n\n-`;
    iter++;
  }while(error > eps && iter < maxIter);
  iterations = `${iterations}\n\nRoot = ${xi}`;
  resultText.value = iterations;
}

//Newton Method
function Newton(xi){
  init();
  let iterations = "";
  xi = parseFloat(xi), iter = 0, error = 100, xiOld = 0;
  do{
    xiOld = xi;
    if(iter != 0){
      xi = r( xi - (f(xi)/fDash(xi)) );
      error = r( Math.abs((xi - xiOld) / xi) * 100 );
    }
    console.log(error);
    iterations = `${iterations}=Iteration: ${iter}=-\nXl: ${xi}\nf(xi): ${f(xi)} | f\`(xi): ${fDash(xi)}\nError%: ${error}%\n\n-`;
    iter++;
  }while (error > eps && iter < maxIter);
  iterations = `${iterations}\n\nRoot = ${xi}`;
  resultText.value = iterations;
}

//Secant Method
function secant(xiMinus1, xi){
  init();
  let iterations = "";
  let error = 0, xiOld = 0, iter = 0;
  xiMinus1 = parseFloat(xiMinus1);
  xi = parseFloat(xi);
  do{
    xiOld = xi;
    if(iter == 0) error = 100;
    else error = r( Math.abs((xi - xiMinus1)/xi) * 100 );
    iterations = `${iterations}-=Iteration: ${iter}=-\nXi-1: ${xiMinus1}\nf(xi-1): ${f(xiMinus1)}\nXi: ${xi}\nf(xi): ${f(xi)}\nError%: ${error}%\n\n`;
    xi = r( xi - ((f(xi) * (xiMinus1 - xi)) / (f(xiMinus1) - f(xi))) );
    console.log(xi);
    xiMinus1 = xiOld;
    iter++;
  }while (error > eps && iter < maxIter)
  iterations = `${iterations}\n\nRoot = ${xiMinus1}`;
  resultText.value = iterations;
}

// ----- CHAPTER 2 -------
var res = ""; // Output string
var partialPivot = false; // Partial Pivot or not?
var swaps = []; // Partiapl Pivot swaps (How many partial pivot swaps?)
var lMatrixSwap = false; // Swap L matrix? (R2 <-> R3)
var step = 1; // Steps counter
var m21, m31, m32;

//Partial Pivot Button Handler
function partial(){
  let partialBtn = document.getElementById("partialBtn");
  if(partialPivot == true){
    partialPivot = false;
    partialBtn.innerHTML = "Partial Pivot = OFF";
  }
  else{
    partialPivot = true;
    partialBtn.innerHTML = "Partial Pivot = ON";
  }
}

function random(){
  for(let i = 0; i < 3; i++){
    for(let j = 0; j < 4; j++){
      let n = Math.floor(Math.random() * 50);
      document.getElementById(`a${i}${j}`).value = n;
    }
  }
}

function setup(id){
  m21 = 0, m31 = 0, m32 = 0, swaps = [], step = 1, lMatrixSwap = false;
  res = "";
  let exam = [[7, 5, 1, 3], [2, 12, 0, 13], [4, 1, 3, 23]];
  let testingPurposes = [[1, 1, -1, 2],[5, 2, 2, 9],[-3, 5, -1, 1]];
  
  let a = [[], [], []];
  
    for(let i = 0; i < 3; i++){
      for(let j = 0; j < 4; j++){
        let v = parseInt(document.getElementById(`m${i}${j}`).value);
        if(isNaN(v)){
          a = exam;
          break;
        }
        else a[i][j] = v;
      }
    }
  if(id == 0) resultText.value = GaussianElimination(a, 0);
  else if(id == 1) resultText.value = LUDecomposition(a);
  else if(id == 3) resultText.value = guassJordan(a);
  else if(id == 2 && partialPivot == false) resultText.value = cramersRule(a);
};

function gaussJordanDivide(matrix, row, divider){
  for(let i = 0; i < 4; i++){
    matrix[row][i] = matrix[row][i]/divider;
  }
  return matrix;
}

function gaussJordanMultiply(matrix, row1, row2){
  let multiplier = -matrix[row1][row2];
  for(let i = 0; i < 4; i++) matrix[row1][i] = matrix[row1][i] + (multiplier * matrix[row2][i]);
  return matrix;
}

function guassJordan(matrix){
  let cols = matrix[0].length;
  let rows = matrix.length;
  // Augmented Matrix

  //Partial pivot or not?
  if(partialPivot == true){
    if(Math.abs(matrix[0][0] < Math.abs(matrix[1][0]) || Math.abs(matrix[0][0]) < Math.abs(matrix[2][0]))){
      if(Math.abs(matrix[1][0]) > Math.abs(matrix[2][0]) || Math.abs(matrix[1][0]) == Math.abs(matrix[2][0])){
        swap(matrix, 0, 1);
      }
      else if(Math.abs(matrix[1][0]) < Math.abs(matrix[2][0])){
        swap(matrix, 0, 2);
      }
    }
  }
  // Make m00 = 1 (divide or swap)
  if(matrix[0][0] != 1){
    if(matrix[1][0] == 1 && partialPivot == false) {
      swap(matrix, 0, 1);
    }
    else if(matrix[2][0] == 1 && partialPivot == false){
      swap(matrix, 0, 2);
    }
    else {
      matrix = gaussJordanDivide(matrix, 0, matrix[0][0]);
    }
  }

  // R2 = R2 + Mutiplier * R1
  if(matrix[1][0] != 0){
    matrix = gaussJordanMultiply(matrix, 1, 0);
  }

  // R3 = R3 + Multiplier * R1
  if(matrix[2][0] != 0){
    matrix = gaussJordanMultiply(matrix, 2, 0);

  // Partial pivot enabled  
  }
  if(partialPivot == true){
    if(Math.abs(matrix[1][1] < Math.abs(matrix[2][1]))){
      swap(matrix, 1, 2);
      lMatrixSwap = true;
    }
  }

  // Make m11 = 1 (divide or swap)
  if(matrix[1][1] != 1){
    if(matrix[2][1] == 1  && partialPivot == false){
      swap(matrix, 1, 2);
    }
    else{
      matrix = gaussJordanDivide(matrix, 1, matrix[1][1]);
    }
  }

  // R3 = R3 + Multiplier R2
  if(matrix[0][1] != 0){
    matrix = gaussJordanMultiply(matrix, 0, 1);
  }
  
  // R1 = R1 + Multiplier R2
  if(matrix[2][1] != 0){
    matrix = gaussJordanMultiply(matrix, 2, 1);
  }

  // Make m22 = 1 (By dividing)
  if(matrix[2][2] != 1){
    matrix = gaussJordanDivide(matrix, 2, matrix[2][2]);
  }

  // R1 = R1 + Multiplier R2
  if(matrix[0][2] != 0){
    matrix = gaussJordanMultiply(matrix, 0, 2);
  }

  // R3 = R3 + Multiplier R3
  if(matrix[2][1] != 0){
    matrix = gaussJordanMultiply(matrix, 1, 2);
  }

  let x3 = matrix[2][3] / matrix[2][2];
  let x2 = (matrix[1][3] - ((matrix[1][2] * x3)))/matrix[1][1];
  let x1 = (matrix[0][3] - ( (matrix[0][1] * x2)  + (matrix[0][2] * x3) ))/matrix[0][0];
  res += `Matrix after elimination:\n${display(matrix)}`;
  res += `\n\nx1: ${x1}\nx2: ${x2}\nx3: ${x3}`;
  return res;
}

function roundTo(num, places) {
  const factor = 10 ** places;
  return Math.round(num * factor) / factor;
};
function reduce(matrix){
  let newMatrix = [[], [], []];
  let cols = matrix[0].length;
  let rows = matrix.length;
  if(cols == 3) return 0;
  for(i = 0; i < cols-1; i++)
    for(j = 0; j < rows; j++)
      newMatrix[i][j] = matrix[i][j]
  return newMatrix;
};

//Displaying the matrix
function display(arr){
  let m = "";
  let cols = arr[0].length;
  if(cols == undefined){
    return `|${arr[0]}|\n|${arr[1]}|\n|${arr[2]}|`;
  }
  let rows = arr.length;

  for(i = 0; i < rows; i++){
    m += `|`
    for(j = 0; j < cols; j++){
      m += ` ${roundTo(arr[i][j], 1)}\t`;
    }
    if(i != rows-1) m += "|\n";
    else m += "|";
  }
  return m;
}

function swap(matrix, r1, r2){
  let cols = matrix[0].length;
  if(cols == undefined){ // While matrix is only one column
    let tempValue = matrix[r1];
    matrix[r1] = matrix[r2];
    matrix[r2] = tempValue;
    return matrix;
  }
  let rows = matrix.length;
  let tempRow = [];
  for(let i = 0; i < cols; i++){
    tempRow[i] = matrix[r1][i];
    matrix[r1][i] = matrix[r2][i];
    matrix[r2][i] = tempRow[i];
  }
  swaps.push([r1, r2]);
  return matrix;
}

// Gaussian Elimination
function GaussianElimination(matrix, op){
  let cols = matrix[0].length;
  let rows = matrix.length;

    //Check for swaps
  if(partialPivot == true){
    if(Math.abs(matrix[0][0] < Math.abs(matrix[1][0]) || Math.abs(matrix[0][0]) < Math.abs(matrix[2][0]))){
      if(Math.abs(matrix[1][0]) > Math.abs(matrix[2][0]) || Math.abs(matrix[1][0]) == Math.abs(matrix[2][0])){
        swap(matrix, 0, 1);
      }
      else if(Math.abs(matrix[1][0]) < Math.abs(matrix[2][0])){
        swap(matrix, 0, 2);
      }
    }
  }
  //Get m21, m31
  m21 = matrix[1][0]/matrix[0][0];
  m31 = matrix[2][0]/matrix[0][0];

  //Applying R2 = R2 - R1(m21)
  for(i = 0; i < cols; i++){
    matrix[1][i] = matrix[1][i] - (matrix[0][i] * m21);
  }

  //Applying R3 = R3 - R1(m31)
  for(i = 0; i < cols; i++){
    matrix[2][i] = matrix[2][i] - (matrix[0][i] * m31);
  }

  //Check for swaps
  if(partialPivot == true){
    if(Math.abs(matrix[1][1] < Math.abs(matrix[2][1]))){
      swap(matrix, 1, 2);
      lMatrixSwap = true;
      }
    }

  //Get m32
  m32 = matrix[2][1]/matrix[1][1];

  //Applying R3 = R3 - R2(m32)
  for(i = 0; i < cols; i++){
    matrix[2][i] = matrix[2][i] - (matrix[1][i] * m32);
  }

  if (op == 1) return matrix;
  else {
    let x3 = matrix[2][3] / matrix[2][2];
    let x2 = (matrix[1][3] - ((matrix[1][2] * x3)))/matrix[1][1];
    let x1 = (matrix[0][3] - ( (matrix[0][1] * x2)  + (matrix[0][2] * x3) ))/matrix[0][0];
    res += `m21: ${m21}\nm31: ${m31}\nm32: ${m32}\n\n`
    res += `Matrix after elimination:\n ${display(matrix)}\n\n`
    res += `x1: ${x1}\nx2: ${x2}\nx3: ${x3}`
    return res;
  }
}

// LU Decomposition
function LUDecomposition(matrix){
  let cols = matrix[0].length;
  let rows = matrix.length;

  //Separate A from B
  let aMatrix = reduce(matrix);
  let bMatrix = [matrix[0][3], matrix[1][3], matrix[2][3]];

  aMatrix = GaussianElimination(aMatrix, 1);

  if(partialPivot == true){
    let num = swaps.length; // Save amount of swaps (Because it will be increased after B and L swaps)
    for(let i = 0; i < num; i++){
      bMatrix = swap(bMatrix, swaps[i][0], swaps[i][1]);
    }
  }
  let c1, c2, c3;
  let lMatrix = [[1, 0, 0],
                [m21, 1, 0],
                [m31, m32, 1]];

  if(lMatrixSwap == true){
    //lMatrix = swap(lMatrix, 1, 2);

    //BUILD
    c1 = bMatrix[0];
    c2 = bMatrix[1] - (lMatrix[1][0] * c1);
    c3 = bMatrix[2] - ((lMatrix[2][0] * c1) + (lMatrix[2][1] * c2));
  } else{
    c1 = bMatrix[0];
    c2 = bMatrix[1] - (lMatrix[1][0] * c1)
    c3 = bMatrix[2] - ((lMatrix[2][0] * c1) + (lMatrix[2][1] * c2));
  
  }
  let uMatrix = matrix;


  // UX = C
  let x3 = c3/aMatrix[2][2];
  let x2 = (c2 - ((aMatrix[1][2] * x3)))/aMatrix[1][1];
  let x1 = (c1 - ( (aMatrix[0][1] * x2)  + (aMatrix[0][2] * x3) ))/aMatrix[0][0];
  res += `m21: ${m21}\nm31: ${m31}\n m32: ${m32}.\n\n`;
  res += `L Matrix:\n${display(lMatrix)}\n\n`;
  res += `c1: ${c1}\nc2 ${c2}\n c3: ${c3}\n\n`;
  res += `U Matrix:\n${display(aMatrix)}\n\n`;
  res += `x1: ${x1}\nx2 ${x2}\n x3: ${x3}\n\n`;
  return res;
}

// All Cramer's Rule Functions
function determine(matrix){
  let D = ( matrix[0][0] * ((matrix[1][1] * matrix[2][2]) - (matrix[1][2] * matrix[2][1])))
        - ( matrix[0][1] * ((matrix[1][0] * matrix[2][2]) - (matrix[1][2] * matrix[2][0])))
        + ( matrix[0][2] * ((matrix[1][0] * matrix[2][1]) - (matrix[1][1] * matrix[2][0])))
           return D;
}
function copyMatrix(a, b){
  for(let i = 0; i < 3; i++)
    for(let j = 0; j < 3; j++)
      a[i][j] = b[i][j]
  return a;
}
function cramersRule(matrix){
  let D, D1, D2, D3;

  //Separate A from B
  let aMatrix = reduce(matrix);
  let bMatrix = [matrix[0][3], matrix[1][3], matrix[2][3]];

  D = determine(aMatrix);

  let tempMatrix = [[], [], []];
  for(let i = 0; i < 3; i++){
    copyMatrix(tempMatrix, aMatrix);
    for(let j = 0; j < 3; j++){
      tempMatrix[j][i] = bMatrix[j]
    }
    if (i == 0){
      D1 = determine(tempMatrix);
      res += `Matrix 1:\n${display(tempMatrix)}\n\n`
    }
    else if (i == 1){
      D2 = determine(tempMatrix);
      res += `Matrix 2:\n${display(tempMatrix)}\n\n`
    }
    else if (i == 2){
      D3 = determine(tempMatrix);
      res += `Matrix 3:\n${display(tempMatrix)}\n\n`
    }
  }

  let x1 = D1/D;
  let x2 = D2/D;
  let x3 = D3/D;
  res += `D: ${D}\nD1: ${D1}\nD2: ${D2}\nD3: ${D3}\n\n`
  res += `x1 = ${D1}/${D} = ${x1}\nx2 = ${D2}/${D} = ${x2}\nx3 = ${D3}/${D} = ${x3}`
  return res;
}

// ----- CHAPTER 3 -------

// Repeat.
function goldenSection(xl, xu, maxIter, goal){
  if (document.getElementById("inputMethod").value != "") fx = document.getElementById("inputMethod").value;
  const R = 0.618;
  let x1, x2, d, iter = 0;
  xl = parseFloat(xl);
  xu = parseFloat(xu)
  maxIter = parseInt(maxIter);

  // Display iterations
  let iterations = "";
  do {
    d = r( R * (xu-xl) );
    x1 = r( xl + d );
    x2 = r( xu - d );
    iterations = `${iterations}=Iteration: ${iter}=-\nxl: ${xl}\nf(xl): ${f(xl)} | x2: ${x2}\nf(x2): ${f(x2)}\nx1: ${x1}\nf(x1): ${f(x1)}\nxu: ${xu}\nf(xu): ${f(xu)}\nd: ${d}\n\n-`;
    if(goal == 'max'){
      if( f(x2) > f(x1) ) xu = x1;
      else if( f(x2) < f(x1) ) xl = x2;
    }else if(goal == 'min'){
      if( f(x2) > f(x1) ) xl = x2;
      else if( f(x2) < f(x1) ) xu = x1;
    }
    iter++;
  } while(iter < maxIter);
  iterations = `${iterations}\n\nD = ${d}`;
  resultText.value = iterations;
}

function gradiant(x, y, expression){
  // Functions
  
  //Define
  let step = 0;
  let output = "";
  let dy, dy2, dydx, dx, dx2, dxdy
  let gradiant;
  let det;
  let result;
  x = parseFloat(x), y = parseFloat(y);
  expression = expression || "2x * y + 2 * x - x^2 - 2 * y^2";
  console.log(expression);
  // Calculate dy and dx
  dx = math.derivative(expression, 'x');
  dy = math.derivative(expression, 'y');
  output += `Step ${step++}: Get df/dx & df/dy\ndf/dx: ${dx.toString()}\ndf/dy: ${dy.toString()}`;

  // Calculate gradiant
  gradiant = r( math.atan(dy.eval({x: x, y: y})/dx.eval({x: x, y: y})) / (Math.PI/180) );
  output += `\n\nStep ${step++}: Calculate gradiant f = df/dx i + df/dy j = tan^-1(${dy.eval({x: x, y: y})}/${dx.eval({x: x, y: y})}) = ${gradiant}`;

  // Output H Matrix
  output += `\n\nStep ${step++}: H matrix:\n| d2f/dx2  d2f/dxdy |\n| d2f/dydx  d2f/dy2 |`;
  dx2 = math.derivative(dx, 'x');
  dxdy = math.derivative(dx, 'y');
  dy2 = math.derivative(dy, 'y');
  dydx = math.derivative(dy, 'x');
  output += `\n\nStep ${step++}:Substitute in H matrix:\n|\t${dx2.eval({x: x, y: y})}\t${dxdy.eval({x: x, y: y})}\t|\n|\t${dydx.eval({x: x, y: y})}\t${dy2.eval({x: x, y: y})}\t|`;

  // Output det(h)
  det = dx2.eval({x: x, y: y}) * dy2.eval({x: x, y: y}) - dxdy.eval({x: x, y: y}) * dydx.eval({x: x, y: y});
  output += `\n\nStep ${step++}: Calculate Det(H)\n| H | = ${dx2.eval({x: x, y: y})} (${dy2.eval({x: x, y: y})}) - ${dxdy.eval({x: x, y: y})} (${dydx.eval({x: x, y: y})}) = ${det}`;

  // Check if local minimum, local maximum, suddle point
  if(det > 0 && dx2 > 0) result = "local max";
  else if(det > 0 && dx2 < 0) result = "local min";
  else if (det <= 0) result = "suddle point";

  output += `\n\nStep ${step++}: Final result: ${result}`;
  resultText.value = output;
}
