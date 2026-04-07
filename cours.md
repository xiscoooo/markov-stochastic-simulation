
# Document de Synthèse : Processus Stochastiques, Files d'Attente et Modèles de Migration

Ce document rassemble de manière exhaustive les concepts liés aux processus stochastiques à temps continu, la théorie des files d'attente, l'épidémiologie (processus de branchement) et les processus de migration fermés (colonies).

---

## 1. Processus de Poisson et Temps Aléatoires

Le processus de Poisson est la brique fondamentale pour modéliser des événements qui se produisent aléatoirement dans le temps (comme des arrivées de clients).

### 1.1. Modélisation Temporelle
* **Instants d'arrivée :** La suite $(T_n)_{n \ge 1}$ représente les instants successifs d'arrivée[cite: 6]. 
* **Temps inter-arrivées :** Le temps écoulé entre deux événements successifs est noté $S_n$[cite: 8]. La relation entre ces variables est $T_n = \sum_{i=1}^n S_i$[cite: 9].
* Les variables $S_n$ sont indépendantes et suivent toutes une loi exponentielle de paramètre $\lambda$, notée $S_n \sim \text{Exp}(\lambda)$[cite: 15, 18].
* L'espérance du temps inter-arrivées est $E(S) = \frac{1}{\lambda}$[cite: 21].
* **Absence de mémoire :** La loi exponentielle vérifie la propriété fondamentale d'absence de mémoire : $P(S > t+s | S > t) = P(S > s) = e^{-\lambda s}$[cite: 22, 24].
* **Loi de l'instant de la n-ième arrivée ($T_n$) :** $T_n$ suit une loi d'Erlang. Sa densité de probabilité est donnée par $f_{T_n}(t) = \frac{\lambda^n}{(n-1)!}e^{-\lambda t}t^{n-1}1_{\mathbb{R}^+}(t)$[cite: 25].
* La transformée de Laplace associée à $T_n$ se calcule via $E(e^{-sT_n}) = \prod_{i=1}^n E(e^{-sS_i}) = \left(\frac{\lambda}{\lambda+s}\right)^n$[cite: 26].

### 1.2. Processus de Comptage $N(t)$
* Le processus $N(t)$ compte le nombre d'événements survenus jusqu'à l'instant $t$[cite: 10, 12].
* $N(t)$ suit une loi de Poisson de paramètre $\lambda t$, c'est-à-dire $P(N(t) = k) = e^{-\lambda t} \frac{(\lambda t)^k}{k!}$[cite: 26].
* Le nombre moyen d'arrivées jusqu'à $t$ est $E(N(t)) = \lambda t$[cite: 27].
* **Propriétés clés :**
    * Les accroissements temporels du processus, par exemple $N(t_1)$ et $N(t_2) - N(t_1)$, sont strictement indépendants[cite: 35].
    * Conditionnellement à l'événement $\{N(t) = n\}$, le vecteur des instants d'arrivée $(T_1, \dots, T_n)$ a la même loi que les statistiques d'ordre $U_{(1)}, \dots, U_{(n)}$ d'un échantillon de $n$ variables uniformes sur l'intervalle $[0, t]$[cite: 37, 40, 44].

### 1.3. Opérations sur les Processus et Paradoxes
* **Superposition :** La somme de deux processus de Poisson indépendants de paramètres respectifs $\lambda_1$ et $\lambda_2$ forme un nouveau processus de Poisson de paramètre $\lambda_1 + \lambda_2$[cite: 54, 56]. Cela s'appuie sur la propriété du minimum : $\min(\text{Exp}(\lambda_1), \text{Exp}(\lambda_2)) \sim \text{Exp}(\lambda_1+\lambda_2)$[cite: 57, 58].
* **Amincissement (Thinning) :** Si l'on filtre chaque événement d'un processus de taux $\lambda$ en le gardant avec une probabilité $p$ et en l'ignorant avec $1-p$, on obtient deux processus de Poisson indépendants de paramètres $\lambda p$ et $\lambda(1-p)$[cite: 59, 61, 64].
* **Paradoxe de l'autobus (Inspection Paradox) :** Si un observateur arrive à un instant $t$ fixe :
    * Le temps d'attente jusqu'au prochain événement, $W_t = T_{N(t)+1} - t$, suit une loi $\text{Exp}(\lambda)$ exacte, comme s'il venait de démarrer[cite: 65, 68, 70]. 
    * Le temps écoulé depuis le dernier événement, $Z_t = t - T_{N(t)}$, suit une loi tronquée $\min(t, \text{Exp}(\lambda))$[cite: 69, 70].
    * Ces deux variables, $W_t$ et $Z_t$, sont indépendantes[cite: 70].

---

## 2. Processus de Markov à Temps Continu (CTMC)


Les processus de Markov à temps continu (CTMC), notés $(X_t, t \ge 0)$, modélisent des systèmes évoluant dans un espace d'états $E$ fini ou dénombrable (souvent $\mathbb{N}$)[cite: 90, 92].

### 2.1. Trajectoires et Générateur Infinitésimal
* **Dynamique de séjour :** En arrivant dans un état $x$, le processus y reste bloqué pendant une durée aléatoire suivant une loi $\text{Exp}(\lambda_x)$[cite: 107, 111]. Il saute ensuite vers l'état $y$ avec une probabilité $q_{x,y}$[cite: 108, 116].
* **Générateur infinitésimal $A$ :** Cette matrice encode complètement la dynamique du système[cite: 216, 217, 221].
    * Les termes extra-diagonaux sont les taux de transition : $a_{x,y} = \lambda_x q_{x,y}$[cite: 218].
    * Le terme diagonal $a_{x,x}$ est strictement négatif tel que la somme de chaque ligne soit nulle : $\sum_{y \neq x} a_{x,y} = -a_{x,x}$[cite: 219, 220].
* **Opérateur de transition :** $A$ opère sur l'espace des fonctions bornées $l^\infty(E) \to l^\infty(E)$[cite: 235]. Pour une fonction $f$, on a $Af(x) = \sum_{y \ne x} (f(y) - f(x)) a_{x,y}$[cite: 248].
* **Formule de Dynkin (Martingale) :** L'espérance de l'évolution d'une fonction $f$ s'écrit : $E(f(X_t)) = E(f(X_0)) + \int_0^t E(Af(X_s)) ds$[cite: 251].

### 2.2. Comportement Transitoire et Stationnaire
* **Régime transitoire :** La distribution de probabilité à un instant $t$ se calcule via l'exponentielle de matrice : $\Pi_t = \Pi_0 e^{tA} = \Pi_0 \sum_{n=0}^{\infty} A^n \frac{t^n}{n!}$[cite: 289, 290].
* **Classification des états :** * Soit $\tau_x = \inf\{t > 0, X_t = x\}$ le temps de premier retour en $x$[cite: 301].
    * Un état $x$ est récurrent si on y retourne presque sûrement ($P_x(\tau_x < +\infty) = 1$)[cite: 302]. Il est transitoire si $P_x(\tau_x = +\infty) > 0$[cite: 303, 304].
    * Un état est récurrent positif si l'espérance du temps de retour est finie ($E_x(\tau_x) < +\infty$)[cite: 308].
    * Si la chaîne est irréductible (graphe fortement connexe), tous les états partagent la même nature[cite: 299, 306].
* **Équilibre Stationnaire :** Pour une chaîne irréductible, l'existence d'une loi stationnaire $\Pi$ (telle que $\Pi A = 0$ et $\sum \Pi(x) = 1$) est équivalente au fait que tous les états soient récurrents positifs[cite: 316, 318]. On a alors $\Pi(x) = \frac{1}{E_x(\tau_x)}$[cite: 319].
* **Théorème Ergodique :** La proportion du temps passé dans un état $y$ sur un intervalle $[0, T]$ converge : $\lim_{T \to \infty} \frac{1}{T} \int_0^T 1_{\{y\}}(X_s) ds = \pi(y)$[cite: 350, 353].

---

## 3. Théorie des Files d'Attente

Le formalisme de Kendall **A/S/s/buf/disc** permet de standardiser la description[cite: 144]:
* **A :** Inter-arrivées (M pour sans mémoire/Exponentielle, D pour Déterministe, GI pour général)[cite: 152, 155, 157].
* **S :** Temps de service[cite: 153].
* **s :** Nombre de serveurs[cite: 145].
* **buf :** Capacité totale (serveurs + file)[cite: 146].
* **disc :** Discipline (FIFO par défaut, LIFO, PS/Processor Sharing)[cite: 147, 150, 151].

### 3.1. Le modèle M/M/1
* Processus d'arrivées $PP(\lambda)$ et temps de service $\text{Exp}(\mu)$[cite: 167, 168].
* **Stabilité :** On pose l'intensité du trafic $\rho = \frac{\lambda}{\mu}$[cite: 335]. Le système est stable et ne sature pas si $\rho < 1$, signifiant que le nombre moyen d'arrivées est strictement inférieur au taux de service possible[cite: 340].
* **Loi stationnaire :** Résolution par équation de balance $\lambda \pi(i) = \mu \pi(i+1)$[cite: 332]. On obtient $\pi(i) = \rho^i (1-\rho)$[cite: 335].
* **Espérance de clients ($Q$) :** Le nombre moyen de clients en régime stationnaire se calcule analytiquement via une dérivation de série : $E(X_\infty) = \sum_{k=0}^{\infty} k(1-\rho)\rho^k = \rho(1-\rho) \frac{d}{d\rho} (\frac{1}{1-\rho}) = \frac{\rho}{1-\rho}$[cite: 363, 365, 369, 388].
* **Loi de Little et Temps d'attente :**
    * La formule de Little affirme $Q = \lambda W$[cite: 382, 385].
    * Le temps de séjour moyen $S$ (attente $W$ + service) est $S = \frac{1}{\mu}\frac{1}{1-\rho}$[cite: 379, 390].
    * Le temps d'attente seul est $W = \frac{1}{\mu}\frac{\rho}{1-\rho}$[cite: 380, 397]. W et S sont robustes face à la discipline de service utilisée (FIFO, LIFO)[cite: 398].

### 3.2. Le modèle M/M/$\infty$ et M/M/s/s (Formule d'Erlang)
* **M/M/$\infty$ :** Une infinité de serveurs fonctionnant à vitesse constante $\mu$. Si $i$ clients sont présents, le taux de service global est $i\mu$. Les taux du générateur sont $\lambda_i = \lambda$ (arrivées) et le taux de départ est $i\mu$[cite: 198, 200, 202].
* **Système à pertes (M/M/s/s) :** Utilisé par exemple pour les réseaux téléphoniques sans attente (appels rejetés si les $s$ serveurs sont occupés)[cite: 430, 447, 451].
* Équations de balance : $\lambda \pi_i = (i+1)\mu \pi_{i+1}$[cite: 442, 443]. On déduit $\pi_i = \frac{\rho^i}{i!} \pi_0$[cite: 444].
* La probabilité de l'état 0 s'obtient par normalisation : $\pi_0 = \left(\sum_{k=0}^s \frac{\rho^k}{k!}\right)^{-1}$[cite: 445].
* **Erlang B :** Le pourcentage d'appels perdus équivaut à la probabilité de l'état bloquant $s$ : $E(\rho, s) = \frac{\rho^s/s!}{\sum_{j=0}^s \rho^j/j!}$[cite: 457, 482].
* Le théorème **PASTA** (Poisson Arrivals See Time Averages) prouve que la proportion empirique des arrivées trouvant le système plein est rigoureusement égale à la moyenne temporelle (la probabilité stationnaire d'être plein)[cite: 479, 480].

---

## 4. Épidémiologie et Propagation de Réseaux



### 4.1. Processus de Branchement (Galton-Watson)
Permet de suivre les arbres de contamination à partir d'un nombre $X_0$ de patients infectés initialement [cite: 558-560].
* Soit $Z_n$ le nombre de personnes infectées à la $n$-ième génération, avec $Z_0 = X_0$[cite: 566, 569]. Soit $X_{n,k}$ le nombre d'individus contaminés par le $k$-ième patient[cite: 568].
* Le taux de reproduction $m = E(Z_1)$ est l'espérance du nombre de contaminations directes causées par une personne[cite: 573].
* **Fonctions génératrices :** $G_{Z_n}(s) = E(s^{Z_n}) = \sum_{k=0}^{\infty} s^k P(Z_n=k)$[cite: 572, 595]. On démontre la composition en cascade : $G_{Z_n} = G_{Z_{n-1}} \circ G_{Z_0} = G_{Z_0} \circ \dots \circ G_{Z_0}$[cite: 574, 589].
* **Probabilité d'extinction :** La probabilité que l'épidémie s'arrête exactement à l'étape $n$ est $u_n = P(Z_n = 0) = G_{Z_n}(0)$[cite: 590, 591].
* $G_{Z_0}$ étant une fonction convexe sur $[0, 1]$[cite: 602, 603]:
    * Si $m \le 1$, l'extinction est certaine ($P(Z_\infty=0) = 1$)[cite: 611, 613, 615].
    * Si $m > 1$, la probabilité d'extinction est une racine strictement inférieure à 1 de l'équation $G_{Z_0}(s) = s$[cite: 616, 617].
* Sur un horizon long, le processus normalisé forme une martingale $Z_n \sim W m^n$ qui converge vers une variable aléatoire non-dégénérée[cite: 621, 623, 625].

### 4.2. Modèles Compartimentaux SIS / SIR
Historiquement formulés en 1937, ils partitionnent la population $M$ en compartiments : Susceptibles ($S$), Infectés ($I$) et Recovered/Guéris ($R$) [cite: 627-630].
* Pour le modèle SIS sans immunité, la conservation de la population impose $S(t) + I(t) = M$[cite: 632, 638]. L'espace des états est $E = \{(s,i) \in \mathbb{N}^2, s+i=M\}$[cite: 639].
* **Taux d'évolution :** Un infecté contacte des individus selon un processus de Poisson. La probabilité de rencontrer un susceptible est $\frac{M-k}{M-1}$, et la probabilité que cette rencontre déclenche l'infection est $\beta$[cite: 641, 644, 645]. La guérison naturelle s'opère à un taux fixe $\alpha$[cite: 651, 664].
* **Politique sanitaire :** Pour éviter que le système ne devienne instable, il faut impérativement $\frac{\beta}{\alpha} < 1$[cite: 670]. On peut agir en diminuant $\beta$ (confinement, gestes barrières) ou en augmentant $\alpha$ (traitements)[cite: 671, 673, 674].
* **Limite macroscopique :** Lorsque la population devient très grande ($M \to \infty$), la variable fractionnaire $X_M(t) = \frac{I_M(t)}{M}$ perd ses fluctuations stochastiques[cite: 679]. La formulation intégrale via le problème de martingale montre que l'évolution de la proportion d'infectés devient asymptotiquement régie par une équation différentielle déterministe[cite: 681, 683].

---

## 5. Processus de Migration Fermés (Colonies)


Les réseaux fermés modélisent un système où un nombre strictement constant $N$ d'individus circule entre $J$ colonies ou stations, sans aucune entrée ou sortie externe[cite: 730, 738].

### 5.1. Formalisme du Réseau
* **Espace des états :** Le système est décrit par un vecteur $n = (n_1, \dots, n_J)$ avec la contrainte de population fixe $S = \{n \in \mathbb{Z}_+^J : \sum_{j=1}^J n_j = N \}$[cite: 738].
* **Opérateur de déplacement $T_{jk}$ :** Un saut correspondant au départ d'un individu de la colonie $j$ vers la colonie $k$ modifie l'état par $T_{jk}n = (n_1, \dots, n_j-1, \dots, n_k+1, \dots, n_J)$[cite: 740, 741]. Ces mouvements surviennent individuellement[cite: 745].
* **Taux factorisés :** Le taux d'une telle transition est assumé séparable : $q(n, T_{jk}n) = \lambda_{jk} \phi_j(n_j)$[cite: 750].
    * $\lambda_{jk}$ représente l'intensité de routage ou la topologie du réseau.
    * $\phi_j(n_j)$ modélise la capacité de service de la file source $j$ (ex: $\min(n, s)$ pour $s$ serveurs, ou $n$ pour des serveurs infinis/migration linéaire)[cite: 751, 754, 755].

### 5.2. Résolution et Équilibre
* **Équations de Trafic :** En imaginant un client unique ($N=1$) évoluant dans le réseau, les flux d'équilibres globaux par station $\alpha_j$ satisfont les équations nodales : $\alpha_j \sum_k \lambda_{jk} = \sum_k \alpha_k \lambda_{kj}$[cite: 759, 762].
* **Structure en forme produit :** Par le Théorème 2.4, la loi stationnaire d'un processus de migration fermé s'écrit de façon remarquablement factorisée : 
  $$\pi(n) = \frac{1}{G_N} \prod_{j=1}^J \frac{\alpha_j^{n_j}}{\prod_{r=1}^{n_j} \phi_j(r)}$$
  où $G_N$ est la constante de normalisation dépendante de la taille de population[cite: 766, 767, 770].
* Cette forme produit se valide via les **équations de balance partielle**. Au lieu de juste compenser le flux total sortant de l'état, ces équations montrent que le flux sortant attribuable à la station $j$ équivaut exactement au flux entrant ciblant spécifiquement cette même station $j$ [cite: 774-776, 788].

### 5.3. Cas Pratique : Le Restaurant
* **Contexte :** Un restaurant compte $N$ tables occupées en permanence (processus fermé). Les convives passent par 3 colonies cycliques : le maître d'hôtel (commande plat, durée $\sim \text{Exp}(\mu_1)$), le sommelier (commande vin, durée $\sim \text{Exp}(\mu_2)$), et le repas libre (durée $\sim \text{Exp}(\nu)$)[cite: 691, 693, 694].
* L'espace des états est $E_N$, où $n_1 + n_2 + n_3 = N$[cite: 696, 697, 699, 702].
* Les phases 1 et 2 ont un seul serveur, avec $\phi_1(n_1) = 1_{\{n_1 > 0\}}$ et $\phi_2(n_2) = 1_{\{n_2 > 0\}}$, tandis que la phase repas a des serveurs infinis ($\phi_3(n_3) = n_3$)[cite: 705, 812].
* La loi stationnaire du réseau est : $\pi(n_1, n_2, n_3) = \frac{1}{G(N)} \alpha_1^{n_1} \alpha_2^{n_2} \frac{\alpha_3^{n_3}}{n_3!}$[cite: 710]. La détermination exacte des constantes de trafic $\alpha_i$ (par exemple en posant conventionnellement $\alpha_1=1$) mène à l'unicité de cette distribution de probabilité[cite: 717, 718].
* On peut ainsi dériver analytiquement que la probabilité conjointe que les deux serveurs soient occupés vaut : $\frac{G(N-2)}{G(N)} \alpha_1 \alpha_2$[cite: 719, 720].