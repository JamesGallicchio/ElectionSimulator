import scala.util.Random

object LinearAlgebra {
  type Scalar = Double

  type Dim = Singleton with Int

  trait Matrix[M <: Dim, N <: Dim] { self =>
    def m: M; def n: N

    def apply(i: Int, j: Int): Scalar

    def row(i: Int): Matrix[1, N] = new Matrix[1, N] { val m: 1 = 1; val n: N = self.n
      override def apply(_i: Int, j: Int): Scalar =
        if (_i == 0) self(i, j) else throw new IllegalArgumentException
    }
    def col(j: Int): Matrix[M, 1] = new Matrix[M, 1] { val m: M = self.m; val n: 1 = 1
      override def apply(i: Int, _j: Int): Scalar =
        if (_j == 0) self(i, j) else throw new IllegalArgumentException
    }
    def elems: Seq[Scalar] = for (i <- 0 until m; j <- 0 until n) yield self(i, j)

    def +(o: Matrix[M, N]): Matrix[M, N] = new Matrix[M, N] { val m: M = self.m; val n: N = self.n
      override def apply(i: Int, j: Int): Scalar = self(i, j) + o(i, j)
    }
    def -(o: Matrix[M, N]): Matrix[M, N] = new Matrix[M, N] { val m: M = self.m; val n: N = self.n
      override def apply(i: Int, j: Int): Scalar = self(i, j) - o(i, j)
    }
    def *[P <: Dim](o: Matrix[N, P]): Matrix[M, P] = new Matrix[M, P] { val m: M = self.m; val n: P = o.n
      override def apply(i: Int, j: Int): Scalar = self.row(i).transpose dot o.col(j)
    }
    def *(s: Scalar): Matrix[M, N] = new Matrix[M, N] { val m: M = self.m; val n: N = self.n
      override def apply(i: Int, j: Int): Scalar = self(i, j) * s
    }

    lazy val transpose: Matrix[N, M] = new Matrix[N, M] { val m: N = self.n; val n: M = self.m
      override def apply(i: Int, j: Int): Scalar = self(j, i)
    }
    lazy val cholesky: Matrix[M, N] = {
      val L = ArrMatrix.empty(m, n)

      for (i <- 0 until m; j <- 0 to i) {
        L(i, j) =
          if (i == j)
            Math.sqrt(self(j,j) - (0 until j).map(k => L(j,k)*L(j,k)).sum)
          else
            1.0/L(j,j) * (self(i,j) - (0 until j).map(k => L(i,k)*L(j,k)).sum)
      }

      L
    }

    lazy val det: Scalar = {
      val L = self.cholesky
      val diag = (0 until Math.min(m, n)).map(i => L(i,i)).product
      diag*diag
    }
    lazy val elemSum: Scalar = elems.sum

    def rows: Seq[Matrix[1, N]] = (0 until m) map row
    def cols: Seq[Matrix[M, 1]] = (0 until n) map col

    override def toString: String = s"Matrix(${
      self.rows.map(_.cols.map(MatrixToScalar).mkString(",")).mkString(";")
    })"
  }
  object Matrix {
    def apply[M <: Dim, N <: Dim](m: M, n: N, elems: Seq[Seq[Scalar]]): Matrix[M,N] = {
      if (elems.size == m && elems.forall(_.size == n))
        new ArrMatrix(m, n, elems.flatten.toArray)
      else throw new IllegalArgumentException
    }

    def zero[M <: Dim, N <: Dim](_m: M, _n: N): Matrix[M, N] = new Matrix[M, N] {
      val m: M = _m; val n: N = _n
      override def apply(i: Int, j: Int): Scalar =
        if (0 <= i && i < m && 0 <= j && j < n) 0.0
        else throw new IllegalArgumentException()
    }

    def identity[N <: Dim](_n: N): Matrix[N, N] = new Matrix[N, N] {
      val m: N = _n; val n: N = _n
      override def apply(i: Int, j: Int): Scalar =
        if (0 <= i && i < m && 0 <= j && j < n)
          if (i == j) 1.0 else 0.0
        else throw new IllegalArgumentException()
    }
  }

  type Vector[N <: Dim] = Matrix[N, 1]
  implicit class VectorOps[N <: Dim](val v: Vector[N]) extends AnyVal {
    def apply(i: Int): Scalar = v.apply(i, 0)
    def dot(o: Vector[N]): Scalar = (0 until v.m).map(i => v(i) * o(i)).sum

    def mag2: Scalar = v dot v
    def normalize: Vector[N] = v * (1.0/mag2)

    def boundElems(low: Scalar, high: Scalar): Vector[N] =
      Vector(v.m, v.elems.map(e => if (e < low) low else if (e > high) high else e))
    def stochastic: Vector[N] = boundElems(0.0, 1.0) * (1.0/v.elemSum)
  }
  object Vector {
    def apply[N <: Dim](n: N, elems: Seq[Scalar]): Vector[N] =
      if (elems.size != n)
        throw new IllegalArgumentException
      else new ArrMatrix(n, 1, elems.toArray)

    def zero[N <: Dim](n: N): Vector[N] = Matrix.zero(n, 1)
  }

  implicit def MatrixToScalar: Matrix[1, 1] => Scalar = _(0,0)

  private class ArrMatrix[M <: Dim, N <: Dim] (val m: M, val n: N,
                                               backing: Array[Double]) extends Matrix[M, N] {
    private[LinearAlgebra] def update(i: Int, j: Int, s: Scalar): Unit = backing(i*n + j) = s
    def apply(i: Int, j: Int): Scalar = backing(i*n + j)
  }
  private object ArrMatrix {
    def empty[M <: Dim, N <: Dim](m: M, n: N): ArrMatrix[M, N] =
      new ArrMatrix[M,N](m, n, new Array[Double](m*n))
  }
}

object Stats {
  import LinearAlgebra._

  trait Sampleable[T] {
    def sample: T
  }

  case class MultiNormal[N <: Dim](n: N, mu: Vector[N], sigma: Matrix[N, N], epsilon: Double = 0.000001) extends Sampleable[Vector[N]] {
    private val L: Matrix[N, N] = (sigma + Matrix.identity[N](n)*epsilon).cholesky
    def sample: Vector[N] = L * MultiNormal.randVec[N](n) + mu
  }
  object MultiNormal {
    def fromSamplesWeighted[N <: Dim](n: N, samples: Seq[(Vector[N], Scalar)]): MultiNormal[N] = {
      val totalWeight: Scalar = samples.map(_._2).sum
      val normSamples = samples.map{ case (v, s) => (v, s/totalWeight)}
      val weightsSquared: Scalar = normSamples.map(_._2).map(Math.pow(_, 2)).sum

      val mu: Vector[N] = normSamples.foldLeft(Vector.zero[N](n)) {
        case (sum, (v, weight)) => sum + (v*weight): Vector[N]
      }

      val sigma: Matrix[N, N] = normSamples.foldLeft(Matrix.zero[N,N](n, n)) {
        case (sum, (v, weight)) =>
          val diff: Vector[N] = v - mu
          sum + diff * diff.transpose * weight
      } * (1.0/(1.0 - weightsSquared))

      new MultiNormal(n, mu, sigma)
    }

    def randVec[N <: Dim](n: N): Vector[N] = Vector(n, (1 to n).map(_ => Random.nextGaussian()))
  }
}
