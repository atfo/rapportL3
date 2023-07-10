La génération de la seconde harmonique a lieu dans un cristal doubleur de fréquence avec des propriétés non-linéaires favorables (figure \ref{fig:cristal}).
La propagation d'ondes électromagnétiques dans un milieu non linéaire donne lieu à une équation d'onde avec un terme source, qui s'écrit (cf annexe \ref{NL}) \ncite{boyd,joffre}
\[ \nabla^2 \v E_q(\v r) + \frac{\omega_q^2}{c^2}\tens\epsilon^{(1)}(\omega_q)\cdot \v E_q(\v r) = - \frac{\omega_q^2}{\epsilon_0 c^2} \v P^\mathsc{NL}_q(\v r) \]
dans le domaine de Fourier en temps avec la convention $\v E(\v r, t) = \mathfrak{Re} \left\{ \sum_{n \in \mathbb Z} \v {\boldsymbol{\mathcal E}}_q (\v r) \e{-i q \omega t} \right\}$,
%, avec n un indice désignant la composante spectrale (il s'agira dans notre cas de l'ordre de l'harmonique), 
avec $\tens \epsilon^{(1)}$ le tenseur de permittivité diélectrique relative associé à la partie linéaire de la polarisation et $\v P^\mathsc{NL}_q$ la partie non-linéaire de la polarisation.

On fait l'hypothèse simplificatrice d'une polarisation linéaire selon un des axes principaux et on travaillera par la suite avec des $\mathcal E_q$ scalaires. En particulier, la permittivité tensorielle $\tens \varepsilon^{(1)}$ est remplacée par sa valeur propre correspondante et donne lieu à un indice optique $n_q^2 = 1 + \varepsilon^{(1)}(\omega)/\varepsilon_0$. En écrivant les différentes harmoniques sous la forme d'une ``onde monochromatique d'amplitude variable'' $\mathcal E_q = \A_q(x,y,z) \e{ik_qz}$ avec $\boxed{k_q =\frac{n_q \omega_q}{c}}$ et en se plaçant dans l'approximation paraxiale, i.e. de lente variation de $\A$ avec $z$ $\left(\frac{\partial^2 \A_q}{\partial z^2} \ll k_q \frac{\partial \A_q}{\partial z}\right)$, on arrive aux équations \ncite{joffre}
\begin{align*}  
	\left\{\v\nabla_\bot + 2 i k_q \frac{\partial}{\partial z} \right\} \A_q =
\end{align*}
avec $\v\nabla_\bot = \pdv{}{x^2} + \pdv{}{y^2}$.

En particulier, pour un champ électrique avec une composante à $\omega$ et une à $2\omega$, le terme quadratique dans la polarisation est donné par \ncite{joffre}
\begin{align*}
	P^{(2)} &= \varepsilon_0 \chi^{(2)} \v E^2  \text{ avec $\chi^{(2)}$ la susceptibilité d'odre 2}\\
	&= \frac{\varepsilon_0 \chi^{(2)}}{4} \left\{\mathcal E_1 e^{-i\omega t} + \mathcal E_1^*i e^{i\omega t} + \mathcal E_2 e^{-2i\omega t} + \mathcal E_2^* e^{2i\omega t} \vphantom{\frac12}\right\}^2 \\
	&= \frac{\varepsilon_0 \chi^{(2)}}{4} \left\{\mathcal E_1^2 \e{-2i\omega t} + \E_1^{*2} \e{2i\omega t} + 2|\mathcal E_1|^2 + 2 \E_1\E_2^* \e{i\omega t} + 2 \E_1^*\E_2 \e{-i\omega t} \vphantom{\frac12}\right. \\
	&\qquad\qquad\qquad \left. + 2 \E_1\E_2 \e{-3 i\omega t}  + 2 \E_1^*\E_2^* \e{3 i\omega t} + \mathcal O(\E_2^2) \vphantom{\frac12}\right\}
\end{align*}
où l'on a considéré que $\chi^{(2)}$ est approximativement le même pour toutes les harmoniques.

On voit donc que l'onde incidente (à $\omega$) conduit à un terme à la fréquence double dans la polarisation, qui va conduire à la création d'une onde à cette seconde harmonique comme souhaité. Cette dernière va conduire à un terme à la fréquence fondamentale qui va affecter l'onde incidente ainsi qu'à un terme à la fréquence triple qui va conduire à une onde à la troisième harmonique et ainsi de suite. 
Dans l'hypothèse où la seconde harmonique est d'amplitude faible par rapport au faisceau incident, nous pouvons cependant négliger ces termes d'ordre supérieur. Cette hypothèse, connue sous le nom \textbf{d'hypothèse de non-déplétion}, est discutée en annexe \ref{ndepl}. Nous négligerons également le terme constant dit de redressement. %Nous nous plaçons cependant dans l'hypothèse que la seconde harmonique 

%En écrivant les différents champs 

%Plus précisément, en écrivant les équations d'onde paraxiales pour les différentes harmoniques, 

Ceci conduit à l'équation d'évolution de l'amplitude de l'onde générée $\A_2$ suivante:
\[ \left\{\v\nabla_\bot + 2 i k_2 \frac{\partial}{\partial z} \right\} \A_2 = - \frac{2 \chi^{(2)} \omega^2}{c^2} \A_1^2 \e{- i (k_2 - 2k_1) z} \]

\subsection{Le problème de l'accord de phase}

Cette équation est beaucoup plus abordable dans l'approximation d'une onde plane avec $\A_1$ et $\A_2$ ne dépendant que de $z$ (et donc $\A_1$ constante dans \textbf{l'hypothèse de non-déplétion}):
\begin{align}
	\dv{\A_2}{z} &= i\frac{\chi^{(2)}\omega}{2cn_2} \A_1^2 \e{-i \Delta k z} \text{ avec } \boxed{\Delta k = k_2 - 2k_1} \label{eq:pwe} \\
	\text{soit } \A_2(L) &= i \frac{\chi^{(2)}\omega}{2 cn_2} \A_1^2 L \operatorname{sinc}\left( \frac{\Delta k L}{2} \right) \e{-i\frac{\Delta k L}{2}}
\end{align}
en sortie du cristal en $z=L$ avec $\A_2(0)=0$ à l'entrée.

Cette équation se comprend très bien en considérant que le carré de l'onde incidente, de nombre d'onde $2k_1$, génère un rayonnement qui se déplace ensuite à $k_2$, de sorte qu'en un point $z$ une onde de phase $2k_1z$ vient s'ajouter à une onde se propageant à $k_2$. Ainsi, l'onde à la seconde harmonique générée en $z$ aura à la sortie du cristal en $z=L$ la phase $2k_1 z + k_2 (L-z) = k_2 L - \Delta k z$ (modulo une constante universelle), de sorte que les ondes générées en $z$ et en $z+L_\mathsc{coh} =z + \frac{\pi}{\Delta k}$ interfèrent destructivement, ce qui conduit à une amplitude nulle en sortie, comme illustré figure \ref{fig:agen}.

Avec le cristal de niobate de lithium utilisé, $L_\mathsc{coh}$ vaut à peine $\SI{3}{\micro\meter}$ alors que le cristal fait $\SI{2}{cm}$.
Afin d'avoir une génération efficace, il faudrait donc $\Delta k = \frac{2\pi}{\lambda_2}(n_2-n_1) = 0$ (condition d'accord de phase) avec $\lambda_2=$\lmbd{532} la longueur d'onde dans le vide de la seconde harmonqiue, soit $n_2 = n_1$. Cela n'est \textit{a priori} pas possible sans dispersion anormale. Une solution consiste alors à exploiter la biréfringence du cristal, mais cette méthode est difficile à réaliser et le $\chi^{(2)}$ correspondant à la polarisation requise est souvent assez faible. 

Nous avons choisi une autre solution qui consiste à fabriquer un cristal dont le $\chi^{(2)}$ varie spatialement. En effet, si $\chi^{(2)}(z) = \chi^{(2)}_0 \e{i k_\chi z}$, cela revient à remplacer $\Delta k$ par $\Delta k_\mathsc{eff} = \Delta k - k_\chi$ dans (\ref{eq:pwe}). En pratique, il est difficile de fabriquer un tel cristal, et on préfère inverser le signe de $\chi^{(2)}$ en inversant l'axe extraordinaire d'un matériau ferroélectrique avec une période $\Lambda$ (cf figure \ref{fig:QPM}). On parle alors de quasi-accord de phase et un tel cristal est dit périodiquement pôlé. Dans ce cas, la décomposition de Fourier de $\chi^{(2)}(z) = \chi^{(2)}_0 \operatorname{sgn}[\cos(2\pi z/ \Lambda)]$ montre que le terme de plus grande amplitude est le fondamental d'amplitude $\chi_\mathsc{eff} = \frac2\pi \chi^{(2)}_0$ et de nombre d'onde $k_\chi = \frac{2\pi}{\Lambda}$. Si $k_\chi = \Delta k$, le terme source correspondant sera accordé en phase sur toute la longueur du cristal et permettra donc de générer la seconde harmonique. Par la suite, on ne tiendra compte que de ce terme, oubliant les autres harmoniques de $\chi^{(2)}(z)$.
Ceci nous conduit à une amplitude 
$$
\begin{align}
	&\A_2(L) = i \frac{\chi_\mathsc{eff} \omega}{2 cn_2} \A_1^2L \\	
	\text{et } i &\frac{\chi_\mathsc{eff}\omega}{2 cn_2} \A_1^2 L \operatorname{sinc}\left( \frac{\Delta k_\mathsc{eff} L}{2} \right) 
	\e{-i\frac{\Delta k_\mathsc{eff} L}{2}}
\end{align}
$$
si le quasi-accord de phase n'est pas respecté.

En termes de puissance, la puissance de la seconde harmonique est quadratique en la puissance du fondamental, avec une efficacité (normalisée)
$$
\begin{align}
	\alpha &= \frac{\mathcal P_2}{\mathcal P_1^2} = \frac{\frac12 n_2 \varepsilon_0 c |\A_2|^2}{\left(\frac12 n_1 \varepsilon_0 c |\A_1|^2\right)^2} 
	= \frac{\chi_\mathsc{eff}^2\omega^2}{2 n_2 n_1^2 \varepsilon_0 c^3} L^2
\end{align}
$$
On notera en particulier la dépendance quadratique en la longueur du cristal, qui serait perdue en l'absence d'accord de phase.

\subsection{Cas des faisceaux gaussiens}

Maintenant que nous avons éclairci l'importance de l'accord de phase, nous pouvons étudier l'effet de l'extension limitée du faisceau. Le faisceau incident produit par le laser est un faisceau gaussien, solution de l'équation de Helmholtz paraxiale \ref{eq:NL} sans terme:
$$
\A_1(r,z) = \sqrt{\frac{2I_{1}}{\pi}} \frac{1}{w_{0}} \frac{1}{1+i\zeta} \e{\frac{-r^{2}}{w_{0}^{2} (1+i\zeta) }}
$$
avec $r=\sqrt{x^2+y^2}$, $w_0$ le waist du faisceau et $\zeta = \frac{z}{z_\mathsc{R}}$ où $z_{R}= n_1 \frac{\pi w_0^2}{\lambda_1}$ une fois dans le cristal d'indice $n_1$ à la fréquence du fondamental.

La seconde harmonique est générée par une onde $\v E_1^2 \propto \left( \e{-\left(\frac{r}{w(z)}\right)^2} \right)^2 = \e{-\left(\frac{r}{w(z)/ \sqrt 2}\right)^2}$, ce qui conduit à la supposition que le faisceau vert produit est ``gaussien'' (le profil transverse est gaussien mais la puissance varie longitudinalement) avec la même longueur de Rayleigh et un waist $\sqrt 2$ fois plus petit (ce qui est cohérent avec la relation entre waist et longueur de Rayleigh pour un faisceau de longuer d'onde moitié, si ce n'est la différence d'indice optique).

On pose donc l'Ansatz suivant $$\A_2(r,z) = \frac{A_2(z)}{1+i\zeta}\e{\frac{-2r^{2}}{w_{0}^{2} (1+i\zeta) }} $$où l'on a toujours $\zeta = (z\lambda_1)/(n_1 \pi w_{0}^{2})$. % Cet Ansatz peut aussi être vu comme une variation de la constante dans la solution de l'équation linéaire (sans ``terme source''). (Pas vrai à cause de n_1 au lieu de n_2!)

$A_2$ vérifie alors l'équation d'évolution suivante

\dots

On obtient donc 
\begin{align}
	\alpha = 	
\end{align}
