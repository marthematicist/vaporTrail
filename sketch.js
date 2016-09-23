var setupGlobalVariables = function() {
  // set canvas size to fill the window
  xRes = windowWidth;
  yRes = windowHeight;
  minRes = min( xRes , yRes );
  maxRes = max( xRes , yRes );
  
  // size of the "field" extends past the edges of the canvas.
  xMin = -0.2*xRes;
  xMax = 1.2*xRes;
  yMin = -0.2*yRes;
  yMax = 1.2*yRes;
  edgeWidth = xMax * 0.1;
  
  
  // Threshold distance. Lines will be draw between particles
  // colser than this value
  distThreshold = 400;
  // Distance at which lines begin to fade
  fadeThreshold = 300;
  // controls how "bright" the image is
  alphaFactor = 0.05;
  
  // R, G, B and alpha values for line color
  lineR = 255;
  lineG = 255;
  lineB = 255;
  lineAlpha = 255;
  lineColor = color( lineR , lineG , lineB , lineAlpha );
  // background color alpha
  bgAlpha = 0;
  
  // time between frames, for physics simulation
  dt = 1.0 / ( 40 );
  
  // constands for physics simulation
  edgeSpringConstant = 50000;
  frictionConstant = 0.01;
  universalConstant = 25000;
  epsilon = 10;

  // number of particles - "dots"
  numDots = 50;
  
  // values for randomizing initial particle velocities
  minVel = 0.05*minRes;
  maxVel = 0.1*minRes;
  // values for randomizing initial particle locaitons
  minDistance = 0.1*minRes;
  
  // average mass of particles
  avgMass = 50;

  
  // 
  frameCounter = 0;
  generateNew = false;
  
  clearFirstTime = true;
  startTime = 0;
  waitTime = 3000;
}

// class definition for Dots
class Dots{
  // constructor
  constructor( N ) {
    this.N = N;
    this.X = new Array( this.N );
    this.V = new Array( this.N );
    this.A = new Array( this.N );
    this.M = new Array( this.N );
    var nSquared = this.N*this.N;
    this.D = new Array(  );
    
    for( var i = 0 ; i < this.N ; i++ ) {
      var farEnough = false;
      while( !farEnough ) {
        this.X[i] = createVector( random(xMin,xMax) , random(yMin , yMax) );
        farEnough = true;
        for( var j = 0 ; j < i ; j++ ) {
          if( p5.Vector.dist( this.X[i] , this.X[j] ) < minDistance ) {
            farEnough = false;
          }
        }
      }
      this.V[i] = p5.Vector.random2D();
      this.V[i].mult( random(minVel,maxVel) );
      this.A[i] = ( createVector( 0 , 0 ) );
      this.M[i] = ( avgMass );
    }
  }
}
    
// method for Dots: updates the distance matrix
Dots.prototype.updateDistances = function() {
  for( var i = 0 ; i < this.N - 1 ; i++ ) {
    for( var j = i ; j < this.N ; j++ ) {
      var d = this.X[i].dist( this.X[j] ); 
      this.D[i*this.N + j] = d; 
      this.D[j*this.N + 1] = d;
    }
  }
}

// method for Dots: retrieves the distance between dots i and j
Dots.prototype.getDist = function( i , j ) {
  var d = this.D[i*this.N + j];
  return d;
}

// method for Dots: zeros out accenerations
Dots.prototype.zeroAccelerations = function() {
  for( var i = 0 ; i < this.N ; i++ ) {
    this.A[i] = createVector( 0 , 0 );
  }
}

//method for Dots: applies friction forces to all dots
Dots.prototype.applyFrictionForces = function() {
  for( var i = 0 ; i < this.N ; i++ ) {
    dA = createVector( this.V[i].x , this.V[i].y );
    dA.mult( -frictionConstant );
    this.A[i].add( dA );
  }
}

// method for Dots: applies repulsivce forces between all dots
Dots.prototype.applyMutualForces = function() {
  for( var i = 0 ; i < this.N - 1 ; i++ ) {
    for( var j = i + 1 ; j < this.N ; j++ ) {
      var d = this.getDist( i , j );
      var f = universalConstant / pow( d * d + epsilon * epsilon , 1.5 );
      var dAi = p5.Vector.sub( this.X[i] , this.X[j] );
      dAi.normalize();
      var dAj = createVector( dAi.x , dAi.y );
      dAi.mult( f * this.M[j] );
      dAj.mult( -f * this.M[i] );
      this.A[i].add( dAi );
      this.A[j].add( dAj ); 
    }
  }
}

// method for Dots: applies edge forces (springy bouncy)
Dots.prototype.applyEdgeForces = function() {
  for( var i = 0 ; i < this.N ; i++ ) {
    var x = this.X[i].x;
    var y = this.X[i].y;
    var m = this.M[i];
    if( x < xMin ) {
        var f = edgeSpringConstant * -x;
        var dA = createVector( f / m , 0 );
        this.A[i].add( dA );
    }
    if( y < yMin ) {
        var f = edgeSpringConstant * ( -y );
        var dA = createVector( 0, f / m );
        this.A[i].add( dA );
    }
    if( x > xMax ) {
        var f = edgeSpringConstant * ( x - ( xMax ) );
        var dA = createVector( -f / m , 0 );
        this.A[i].add( dA );
    }
    if( y > yMax ) {
        var f = edgeSpringConstant * ( y - ( yMax ) );
        var dA = createVector( 0 , -f / m );
        this.A[i].add( dA );
    }
  }
}

// evolves the physics simulation one half step (half of dt)
Dots.prototype.evolveHalfStep = function() {
  this.zeroAccelerations();
  this.updateDistances();
  this.applyMutualForces();
  this.applyEdgeForces();
  this.applyFrictionForces();
  for( var i = 0 ; i < this.N ; i++ ) {
    
    this.V[i].add( p5.Vector.mult( this.A[i] , dt / 2 ) );
  }
}

// method for Dots: evolves the physics simulation n steps
Dots.prototype.evolveFullStep = function( num ) {
  for( var n = 0 ; n < num ; n++ ) {
    for( var i = 0 ; i < this.N ; i++ ) {
      this.X[i].add( p5.Vector.mult( this.V[i] , dt ) );
    }
    this.zeroAccelerations();
    this.updateDistances();
    this.applyMutualForces();
    this.applyEdgeForces();
    this.applyFrictionForces();
    for( var i = 0 ; i < this.N ; i++ ) {
      this.V[i].add( p5.Vector.mult( this.A[i] , dt ) );
    }
  }
}

// method for Dots: draws the particles (currently not used)
Dots.prototype.drawDots = function() {
  for( var i = 0 ; i < this.N ; i++ ) {
    var x = this.X[i].x;
    var y = this.X[i].y;
    ellipse( x , y , 10 , 10 );
    //console.log(x , y);
    
  }
}

// method for Dots: draws lines between particles
Dots.prototype.drawDistances = function(  ) {
  for( var i = 0 ; i < this.N - 1 ; i++ ) {
    for( var j = i + 1 ; j < this.N ; j++ ) {
      var d = this.getDist( i , j );
      if( d < distThreshold ) {
        var x1 = this.X[i].x;
        var y1 = this.X[i].y;
        var x2 = this.X[j].x;
        var y2 = this.X[j].y;
        if( d > distThreshold - fadeThreshold ) {
          var fadeAlpha = lineAlpha * (distThreshold - d) / fadeThreshold;
          lineColor = color( lineR , lineG , lineB , fadeAlpha*alphaFactor );
        } else {
          lineColor = color( lineR , lineG , lineB , lineAlpha*alphaFactor );
        }
        stroke( lineColor );
        line( x1 , y1 , x2 , y2 );          
      }
    }
  }
}



setup = function() {
  //console.log( 'Setting up global variables...' );
  setupGlobalVariables();
  createCanvas( xRes , yRes );
  //console.log( 'setting up Dots object' );
  d = new Dots( numDots );
  //console.log( 'evolving half step');
  d.evolveHalfStep();
  background( 0 , 0 , 0 , 255 );
  
  startTime = millis();
  
  // show the title screen
  textAlign( CENTER );
  fill( 255 , 255 , 255 , 255 );
  textSize( 60 );
  text("VAPOR TRAIL" , 0.5*xRes , 0.5*yRes );
  textSize( 30 );
  text( "-marthematicist-" , 0.5*xRes , 0.5*yRes + 35 );
}

draw = function() {
  // if title screen wait time has not passed, do nothing
  if( millis() - startTime < waitTime ) {
    
    return
  }
  // clear the title screen before the first frame
  if( clearFirstTime ) {
    background( 0 , 0 , 0 , 255 );
    clearFirstTime = false;
  }
  
  // sets fill color
  fill( 255 , 255 , 255 , 255 );
  noStroke();

  //console.log( 'draw function: evolving full step' );
  d.evolveFullStep(2);
  
  //console.log( 'drawing dots' );
  //d.drawDots();
  // console.log( 'drawinf distance lines' );
  d.drawDistances(  );

  frameCounter++;
}

// save an image of the canvas if 's' is typed
function keyTyped() {
  if( key === 's' ) {
    saveCanvas( 'canvas' , 'jpg' );
    console.log("saved");
  }
}
