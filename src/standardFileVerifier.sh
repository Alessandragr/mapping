#!/bin/bash
# -*- coding: utf-8 -*-

# Le code doit recevoir une archive et vérifier si c'est bien un fichier archive (et non un dossier) de type `.sam`.  Renvoyer une erreur si c'est le cas.

# Vérifier si cette archive est vide et renvoyer une erreur si c’est le cas. 

# Vérifier si le fichier contient le nombre correct de colonnes pour un fichier `.sam` (minimum 11 colonnes). Émettre une alerte pour un nombre de colonnes supérieur et une erreur si inférieur en disant que c'est pas un archive du type `.sam`.

# Combien de reads sont mappés ? Compter le nombre de reads en fonction du flag (colonne #2).

# Combien de reads sont mappés par chromosome ? Compter les reads par chromosome. 

# Où les reads sont-ils mappés ? Est-ce que l’alignement est homogène le long de la séquence de référence ? Compter le nombre de reads par position sur le chromosome.

# Comment les reads (et paires de reads) sont-ils mappés ? Compter le nombre de reads pour chaque valeur de flag.

# Avec quelle qualité les reads sont-ils mappés ? Compter les reads pour chaque score de qualité ou par tranche de valeurs de score de mapping.

# En fonction des questions 2 et 4, implémenter des paramètres pour filtrer les données. Par exemple, ne conserver que les reads complètement mappés ou ceux avec un score de qualité de mapping inférieur à 30.

# Prendre les reads, réaliser le mapping avec un génome de référence et aligner les résultats.

# EXTRA : Projet additionnel : Coder pour enregistrer le **CIGAR** et stocker ces informations afin de connaître le nombre de délétions, insertions, etc., dans chaque séquence.

# EXTRA : Projet additionnel : Analyser la qualité du CIGAR. Plus un CIGAR est long, plus la qualité du read est incertaine, tandis que les CIGAR plus courts indiquent un alignement plus sûr avec le génome de référence.

# Split - verifier la quantité de colone

#!/bin/bash

#!/bin/bash
echo "Paramètres: $*"
for fichier in "$@"
do
    # Vérification si c'est un fichier régulier et si l'extension est .sam
    if [[ -f "$fichier" && "$fichier" == *.sam ]]
     then
        echo "Le fichier $fichier est un fichier du type .sam"
        
        # Vérifie si le fichier est vide
        if [[ ! -s "$fichier" ]]
         then
            echo "Erreur : Le fichier '$fichier' est vide."
            exit 1
        else
            echo "Le fichier '$fichier' n'est pas vide."
        fi

        # Verification du nombres de colones
        ligne_compteur=0

        # Vérification des colonnes sur les trois premières lignes non-en-tête
        while IFS= read -r ligne || [[ -n "$ligne" ]]
         do
            # Ignorer les lignes d'en-tête (commençant par '@')
            if [[ $ligne == @* ]]
            then
                continue
            fi
            
            # Compte le nombre de colonnes (séparées par des tabulations)
            nb_colonnes=$(echo "$ligne" | awk -F'\t' '{print NF}')

            # Vérifie le nombre de colonnes
            if (( nb_colonnes < 11 ))
            then
                echo "Erreur : La ligne '$ligne' dans '$fichier' contient seulement $nb_colonnes colonnes. Ce n'est pas un fichier .sam valide."
                exit 1
            fi

            # Incrémente le compteur et s'arrête après les 3 premieres lignes
            ((ligne_compteur++))
            if (( ligne_compteur == 3 ))
             then
                break
            fi
        done < "$fichier"

        echo "Le fichier '$fichier' respecte le nombre minimum de colonnes ."

    else
        echo "Ce n'est pas un fichier du type .sam. Insérez un nouveau fichier."
    fi
done

