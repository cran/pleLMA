#' Dateframe of responses to items from depression, anxiety, and stress scales
#'
#' The dass data are responses from a random sample of 1,000 individuals
#' collected during the period 2017 -- 2019. The data were retrieved July
#' 2020.  The items included here are on scales designed to measure
#' depression (14 items), anxiety (13 items), and stress (15 items).  The
#' 42 items were presented online and in random order.  Respondents were
#' instructed to consider the last week when responding to the items using the
#' following categories:
#' (1) Did not apply to me at all;
#' (2) Applied to me to some degree, or some of the time;
#' (3) Applied to me to a considerable degree, or a good part of the time;
#' (4) Applied to me very much, or most of the time.
#'
#' @format A data frame with 1,000 rows (respondents) and 42 columns (items):
#' \describe{
#'   \item{d1}{I couldn't seem to experience any positive feeling at all.}
#'   \item{d2}{I just couldn't seem to get going}
#'   \item{d3}{I felt that I had nothing to look forward to}
#'   \item{d4}{I felt sad and depressed}
#'   \item{d5}{I felt that I had lost interest in just about everything}
#'   \item{d6}{I felt I wasn't worth much as a person}
#'   \item{d7}{I felt that life wasn't worthwhile}
#'   \item{d8}{I couldn't seem to get any enjoyment out of the things I did}
#'   \item{d9}{I felt down-hearted and blue}
#'   \item{d10}{I was unable to become enthusiastic about anything}
#'   \item{d11}{I felt I was pretty worthless}
#'   \item{d12}{I could see nothing in the future to be hopeful about}
#'   \item{d13}{I felt that life was meaningless}
#'   \item{d14}{I found it difficult to work up the initiative to do things}
#'   \item{a1}{I was aware of dryness of my mouth}
#'   \item{a2}{I experienced breathing difficulty (eg, excessively rapid
#'            breathing, breathlessness in the absence of physical exertion)}
#'   \item{a3}{I had a feeling of shakiness (eg, legs going to give way)}
#'   \item{a4}{I felt that I was using a lot of nervous energy}
#'   \item{a5}{I had a feeling of faintness}
#'   \item{a6}{I perspired noticeably (eg, hands sweaty) in the absence of high
#'              temperatures or physical exertion}
#'   \item{a7}{I felt scared without any good reason}
#'   \item{a8}{I had difficulty in swallowing}
#'   \item{a9}{I was aware of the action of my heart in the absence of physical
#'           exertion (eg, sense of heart rate increase, heart missing a beat)}
#'   \item{a10}{I felt I was close to panic}
#'   \item{a11}{I felt terrified}
#'   \item{a12}{I was worried about situations in which I might panic and make a
#'              fool of myself}
#'   \item{a13}{I experienced trembling (eg, in the hands)}
#'   \item{s1}{I found myself getting upset by quite trivial things}
#'   \item{s2}{I tended to over-react to situations}
#'   \item{s3}{I found it difficult to relax}
#'   \item{s4}{I found myself in situations that made me so anxious I was most
#'             relieved when they ended}
#'   \item{s5}{I found myself getting upset rather easily}
#'   \item{s6}{I found myself getting impatient when I was delayed in any way
#'             (eg, elevators, traffic lights, being kept waiting)}
#'   \item{s7}{I felt that I was rather touchy}
#'   \item{s8}{I found it hard to wind down}
#'   \item{s9}{I found that I was very irritable}
#'   \item{s10}{I found it hard to calm down after something upset me}
#'   \item{s11}{I feared that I would be thrownnoff by some trivial but
#'             unfamiliar task}
#'   \item{s12}{I found it difficult to tolerate interruptions to what I was
#'              doing}
#'   \item{s13}{I was in a state of nervous tension}
#'   \item{s14}{I was intolerant of anything that kept me from getting on with
#'              what I was doing}
#'   \item{s15}{I found myself getting agitated}
#' }
#'
#' @source \url{https://openpsychometrics.org}
#'
#' @examples
#' data(dass)
#' head(dass)
#'
"dass"
