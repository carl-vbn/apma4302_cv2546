## Problem 1

**(a)**
$$Au=b$$
$$A\hat{u}=\hat{b}$$
$$A(u+\Delta u)=b+\Delta b$$
$$Au+A\Delta u = b + \Delta b$$
$$b+A \Delta u = b + \Delta b$$
$$A\Delta u = \Delta b$$
$$\Delta u = A^{-1} \Delta b$$
$$||\Delta u|| = ||A^{-1} \Delta b|| \leq ||A^{-1}|| \cdot ||\Delta b ||$$
$$||u - \hat{u}|| \leq ||A^{-1}|| \cdot || b - \hat{b} ||$$
$$\dfrac{||u - \hat{u}||}{||u||} \leq \dfrac{||A^{-1}|| \cdot || b - \hat{b} ||}{||u||}$$

Since $||b||=||Au|| \leq ||A|| \cdot ||u||$, $||u|| \geq \dfrac{||b||}{||A||}$
$$\dfrac{||u - \hat{u}||}{||u||} \leq \dfrac{||A||\cdot||A^{-1}|| \cdot || b - \hat{b} ||}{||b||}$$
$$\boxed{\dfrac{||u - \hat{u}||}{||u||} \leq \kappa(A)\cdot\dfrac{|| b - \hat{b} ||}{||b||}}$$

**(b)**
$$A^{-1}r = A^{-1}Ae$$
$${e=A^{-1}r}$$
$$||e|| = ||A^{-1}r|| \leq ||A^{-1}||\cdot||r||$$
$$\frac{||e||}{||u||} \leq \frac{||A^{-1}||\cdot||r||}{||u||}$$

Since $||b||=||Au|| \leq ||A||\cdot||u||$ so $\dfrac{1}{||u||}\le \dfrac{||A||}{||b||}$:
$$\frac{||e||}{||u||} \leq ||A^{-1}||\cdot||r||\cdot\frac{||A||}{||b||}$$
$${\frac{||e||}{||u||} \leq ||A||\,||A^{-1}||\cdot\frac{||r||}{||b||}}$$
$${\frac{||e||}{||u||} \leq \kappa(A)\cdot\frac{||r||}{||b||}}$$

From $r=Ae$:
$$||r|| = ||Ae|| \leq ||A||\cdot||e||$$
$$||e|| \geq \frac{||r||}{||A||}$$
$$\frac{||e||}{||u||} \geq \frac{||r||}{||A||}\cdot\frac{1}{||u||}$$

Since $u=A^{-1}b$,  $||u||=||A^{-1}b|| \leq ||A^{-1}||\cdot||b||$ so $\dfrac{1}{||u||}\ge \dfrac{1}{||A^{-1}||\cdot||b||}$:
$$\frac{||e||}{||u||} \geq \frac{||r||}{||A||}\cdot\frac{1}{||A^{-1}||\cdot||b||}$$
$${\frac{||e||}{||u||} \geq \frac{1}{||A||\,||A^{-1}||}\cdot\frac{||r||}{||b||}}$$
$${\frac{||e||}{||u||} \geq \frac{1}{\kappa(A)}\cdot\frac{||r||}{||b||}}$$
Thus,

$$\boxed{\frac{1}{\kappa(A)}\frac{||r||}{||b||}\le \frac{||e||}{||u||}\le \kappa(A)\frac{||r||}{||b||}}$$

**Interpretation**: a smaller $k(A)$ (i.e. well conditioned matrix) means the relationship between error $e$ and residual $r$ is tight. For a higher $k(A)$, the relationship is less specific and a small $r$ can mean a large $e$.

## Problem 2

**(a)**

Let
$(v_j)_i=\sin(j\pi x_i)=\sin\!\left(\frac{j\pi i}{m}\right), \qquad i=1,\dots,m-1.$

$$(Av_j)_i=\frac{1}{h^2}\left(-(v_j)_{i-1}+2(v_j)_i-(v_j)_{i+1}\right)$$
$$(Av_j)_i=\frac{1}{h^2}

(-\sin\!(\frac{j\pi(i-1)}{m})
+2\sin\!(\frac{j\pi i}{m})
-\sin\!(\frac{j\pi(i+1)}{m})).$$

Since $\sin(\alpha+\theta)+\sin(\alpha-\theta)=2\sin(\alpha)\cos(\theta),$ with $\alpha=\frac{j\pi i}{m}$ and $\theta=\frac{j\pi}{m}$:

$$\sin\!\left(\frac{j\pi(i+1)}{m}\right)+\sin\!\left(\frac{j\pi(i-1)}{m}\right)
=2\sin\!\left(\frac{j\pi i}{m}\right)\cos\!\left(\frac{j\pi}{m}\right).$$
$$(Av_j)_i=\frac{1}{h^2}
\left(2-2\cos\!\left(\frac{j\pi}{m}\right)\right)
\sin\!\left(\frac{j\pi i}{m}\right).$$

Thus
$$Av_j=\lambda_j v_j,$$
so vectors $\boxed{(v_j)_i=\sin(j\pi x_i)}$ are eigenvectors of $A$.

**(b)**

From (a) we know that ${\lambda_j=\frac{1}{h^2}\left(2-2\cos\!\left(\frac{j\pi}{m}\right)\right).}$

Using $2-2\cos\theta=4\sin^2(\theta/2)$,

$$\boxed{\lambda_j=\frac{4}{h^2}\sin^2\!\left(\frac{j\pi}{2m}\right), \qquad j=1,\dots,m-1.}$$

**(c)**

Since $A$ is symmetric positive definite,
$$\kappa(A)=\frac{\lambda_{\max}}{\lambda_{\min}}.$$

$$\lambda_{\min}=\lambda_1=\frac{4}{h^2}\sin^2\!\left(\frac{\pi}{2m}\right).$$
$$\lambda_{\max}=\lambda_{m-1}
=\frac{4}{h^2}\sin^2\!\left(\frac{(m-1)\pi}{2m}\right)
=\frac{4}{h^2}\cos^2\!\left(\frac{\pi}{2m}\right).$$

So:
$$\kappa(A)=\frac{\cos^2\left(\frac{\pi}{2m}\right)}
{\sin^2\left(\frac{\pi}{2m}\right)}
=\cot^2\left(\frac{\pi}{2m}\right).$$

As $m\to\infty$, $\frac{\pi}{2m}\to 0$ and $\sin \frac{\pi}{2m}\sim \frac{\pi}{2m} \sim 0$, $\cos \frac{\pi}{2m}\sim 1$:

$$\kappa(A)=\cot^2(\frac{\pi}{2m})\sim \left(\frac{1}{\frac{\pi}{2m}}\right)^2
=\left(\frac{2m}{\pi}\right)^2.$$
Thus,
$$\boxed{\kappa(A)=O(m^2)\quad\text{as }m\to\infty.}$$

## Problem 3

$$u(x)=sin(k\pi x) + c (x - \frac{1}{2})^3$$
$$u'(x)=k\pi \cdot cos(k\pi x) + 3c (x - \frac{1}{2})^2$$
$$u''(x)=-k^2\pi^2 \cdot sin(k\pi x) + 6c (x - \frac{1}{2})$$
$$-u''(x)+\gamma u(x)=f(x)$$
$$f(x)=-\left(-k^2\pi^2\sin(k\pi x)+6c\left(x-\frac12\right)\right)+\gamma\left(\sin(k\pi x)+c\left(x-\frac12\right)^3\right)$$
$$f(x)=k^2\pi^2\sin(k\pi x)-6c\left(x-\frac12\right)+\gamma\sin(k\pi x)+\gamma c\left(x-\frac12\right)^3$$
$$f(x)=\left(k^2\pi^2+\gamma\right)\sin(k\pi x)-6c\left(x-\frac12\right)+\gamma c\left(x-\frac12\right)^3$$  
$$\boxed{f(x)=\left(k^2\pi^2+\gamma\right)\sin(k\pi x)-6c\left(x-\frac12\right)+\gamma c\left(x-\frac12\right)^3}$$

## Problem 4

*Note: in the results below, iteration 0 is not counted.*
**(a)** 1000+. Did not meet required tolerances after 1000 iterations.
**(b)** 102 iterations.
**(c)** 1 iteration. This makes sense because without the cubic factor, the residual is enough for CG to descent in one iteration, because the RHS lies entirely in the span of one eigenvector of A.
**(d)** 1 iteration.
**(e)** 1 iteration.
**(f)**  1 iteration for both 1 and 2 processors.