import scala.math.pow

/*Main MAF EM function
 *  IN: Phi - an initial MAF vector of length number of SNPs
 *      GL - Array of arrays of likelihood triples P( D | g )
 *          (note these are NOT multiplied by P(g | phi)! )
 *  OUT: Phi - ML estimate of MAF's across SNPs
 * 
 *  Note: GL is currently an array of (numSnps) arrays of length (numInds),
*/
def emForMAF(Phi: Array[Double], GL: Array[ Array[ (Double, Double, Double) ]]): Array[Double] = {
  var eps = 1.0
  val tol = 0.0001
  L = Phi.length
  while( eps > tol ){
    val Phi_next = Array.fill(L){0.0}
    Phi_next.indices.foreach( i => GL(i).foreach( l =>
      Phi_next(i) += (1.0/(2.0*GL(i).length))*(( 1.0*l._2*2.0*Phi(i)*(1-Phi(i)) + 2.0*l._3*pow(Phi(i),2.0))/( l._1*pow(1.0-Phi(i),2.0) + l._2*2.0*Phi(i)*(1.0-Phi(i)) + l._3*pow(Phi(i),2.0) ))))
    eps = 0.0
    Phi.indices.foreach(i => eps += pow( Phi(i) - Phi_next(i), 2.0 ))
    Phi = Phi_next
  }
  return Phi
}

/*Helper function to compute Y iteratively
For each site, executes the recursion in 4.2.3. Y(i) is Ynk vector for site i
*/
def compY( GL: Array[ Array[ (Double,Double,Double) ]] ): Array[ Array[ Double ]] = {
  L = GL.length
  GLt = GL.transpose
  n = GLt.length
  M = 2.0*n
  Y = Array.ofDim[Double](L,n+1,M+1)
  // NOTE: this ordering may be suboptimal?
  for( i <- 0 to L-1 ){
    for( k <- 0 to M ){
      for( j <- 0 to n ){ // 0 = 0 people not first person
        if( j == 0 ) 
          Y(i)(j)(k) = 1.0
        else if( k == 0 ) 
          Y(i)(j)(k) = (1.0/( 2.0*j*(2.0*j-1.0))) * ( (2.0*j-k)*(2.0*j-k-1.0)*Y(i)(j-1)(k)*GL(i)(j)._1 )
        else if( k == 1 ) 
          Y(i)(j)(k) = (1.0/( 2.0*j*(2.0*j-1.0))) * ( (2.0*j-k)*(2.0*j-k-1.0)*Y(i)(j-1)(k)*GL(i)(j)._1 + 2.0*k*(2.0*j-k)*Y(i)(j-1)(k-1)*GL(i)(j)._2 )
        else 
          Y(i)(j)(k) = (1.0/( 2.0*j*(2.0*j-1.0))) * ( (2.0*j-k)*(2.0*j-k-1.0)*Y(i)(j-1)(k)*GL(i)(j)._1 + 2.0*k*(2.0*j-k)*Y(i)(j-1)(k-1)*GL(i)(j)._2 + k*(k-1.0)*Y(i)(j-1)(k-2)*GL(i)(j)._2 )
      }
    }
  }

  Yr = Array.ofDim[double](L,M)
  for( l <- 0 to L-1 ) Yr(l) = Y(l)(n)
  return Yr
}

/*Main AFS EM function
 *   IN: Phi - an initial MAF vector of length number of SNPs
 *       GL - Array of arrays of likelihood triples P( D | g )
 *           (note these are NOT multiplied by P(g | phi)! )
 *   OUT: Phi - ML estimate of MAF's across SNPs
 *   Note: GL is currently an array of (numSnps) arrays of length (numInds), which is transposed
*/
def emForAFS( Phik: Array[Double], GL: Array[ Array[ (Double, Double, Double) ]]): Array[Double] = {
  val GLt = GL.transpose
  val tol = 0.0001
  val L = GL.length
  val M = Phik.length
  var eps = 1.0
  Y = compY( GL )
  while( eps > tol ){
    var sums = Array.fill(L){0.0}
    sums.indices.foreach( a => Phik.indices.foreach( p => sums(a) += Phik(p)*Y(a)(p) ))
    val Phik_next = Array.fill(M){0.0}
    Phik_next.indices.foreach( i => Y.foreach( y => Phik_next(i) += (1.0/L)*Phik(i)*y(i)/sums(i) ) )
    eps = 0.0
    Phik.indices.foreach( eps += pow( Phik(i) - Phik_next(i), 2.0 ))
    Phik= Phik_next
  }
  return Phik
}
