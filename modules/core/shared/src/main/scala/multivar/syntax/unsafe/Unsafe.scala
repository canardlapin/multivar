package multivar
package syntax
package unsafe

import multivar.core.MultivarError

/** Escape hatch that turns a typed multivar failure into an exception.
  *
  * Prefer matching on `Either` or composing with `flatMap`. Use this only at a
  * deliberate boundary where failure is a programming error rather than a
  * recoverable data condition.
  */
final class MultivarException(val error: MultivarError)
    extends RuntimeException(error.message)

extension [A](value: Either[MultivarError, A])
  def orThrow: A =
    value match
      case Right(success) => success
      case Left(error)    => throw MultivarException(error)
