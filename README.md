# 🧬 Modélisation Plasma par Algorithme Génétique Binaire

Projet numérique réalisé dans le cadre de la Licence EEA (Électronique, Énergie Électrique, Automatique) à l'Université de Toulouse.

---

## 📋 Description du Projet
Ce projet consiste à développer un **algorithme génétique binaire** en langage **C** couplé à une analyse sous **MATLAB** pour résoudre un problème d'interpolation de données expérimentales de plasma. L'objectif est d'ajuster (fitter) des paramètres physiques complexes (tels que $Y_0$, $B$, $\lambda_0$ et $\delta\lambda_L$) sur un profil de raie spectrale en minimisant l'erreur par rapport à des données bruitées et lisses.

---

## 🛠️ Fonctionnalités Clés
- **Codage et Décodage Binaire** : Conversion de variables continues en chromosomes binaires sur 8 bits.
- **Opérations Génétiques** :
  - **Sélection** par tournoi (duels).
  - **Croisement** (crossover) mono-point.
  - **Mutation** aléatoire des gènes basée sur un taux configurable.
- **Évaluation de l'Erreur** : Calcul de l'écart quadratique par rapport au modèle théorique de Lorentz (`F_lorentzienne`).
- **Visualisation** : Exportation des résultats sous format texte pour traçage et analyse sous MATLAB.

---

## 🗂️ Structure du Dépôt
```text
├── RAMMAL_XU_WEI_Projet_Niveau_A_GENETIQUE_CODE.c   # Code source principal en C
├── RAMMAL_XU_WEI_Projet_Niveau_A_Tracee.m            # Script MATLAB pour le tracé des courbes
├── Data_Genetique.txt                              # Fichier de données expérimentales d'entrée
├── Tracee_Finale.txt                               # Fichier de sortie généré par le programme C
├── RAMMAL_XU_WEI_Projet_Niveau_A_GENETIQUE_RAPPORT.pdf # Rapport complet du projet
└── Ziad RAMMAL Presentation.pdf                    # Support de présentation du projet
```

---

## 🚀 Installation et Utilisation

### 1. Prérequis
- Un compilateur C (GCC, Clang, etc.)
- MATLAB (pour l'affichage des courbes)

### 2. Compilation et Exécution (en C)
Ouvre un terminal dans le dossier du projet et exécute les commandes suivantes :

```bash
# Compiler le programme
gcc RAMMAL_XU_WEI_Projet_Niveau_A_GENETIQUE_CODE.c -o projet_genetique -lm

# Exécuter le programme
./projet_genetique
```

### 3. Visualisation des résultats (sous MATLAB)
1. Lance MATLAB.
2. Ouvre et exécute le script `RAMMAL_XU_WEI_Projet_Niveau_A_Tracee.m` pour afficher la superposition des données expérimentales et de la courbe théorique ajustée.

---

## 👥 Auteurs
Projet développé en groupe par :
- **Ziad RAMMAL**
- **Xing WEI**
- **Jinhao XU**

*Université de Toulouse — Licence EEA*