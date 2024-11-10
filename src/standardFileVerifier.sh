#!/bin/bash
# -*- coding: utf-8 -*-


# Combien de reads sont mappés ? Compter le nombre de reads en fonction du flag (colonne #2).

# Combien de reads sont mappés par chromosome ? Compter les reads par chromosome. 

# Où les reads sont-ils mappés ? Est-ce que l’alignement est homogène le long de la séquence de référence ? Compter le nombre de reads par position sur le chromosome.

# Comment les reads (et paires de reads) sont-ils mappés ? Compter le nombre de reads pour chaque valeur de flag.

# Avec quelle qualité les reads sont-ils mappés ? Compter les reads pour chaque score de qualité ou par tranche de valeurs de score de mapping.

# En fonction des questions 2 et 4, implémenter des paramètres pour filtrer les données. Par exemple, ne conserver que les reads complètement mappés ou ceux avec un score de qualité de mapping inférieur à 30.

# Prendre les reads, réaliser le mapping avec un génome de référence et aligner les résultats.

# EXTRA : Projet additionnel : Coder pour enregistrer le **CIGAR** et stocker ces informations afin de connaître le nombre de délétions, insertions, etc., dans chaque séquence.

# EXTRA : Projet additionnel : Analyser la qualité du CIGAR. Plus un CIGAR est long, plus la qualité du read est incertaine, tandis que les CIGAR plus courts indiquent un alignement plus sûr avec le génome de référence.

echo "Paramètres: $*"
for fichier in "$@"
do
    # Vérification si c'est un fichier régulier et si l'extension est .sam
    if [[ -f "$fichier" && "$fichier" == *.sam && -r "$fichier" ]]
     then
        echo "Le fichier $fichier est un fichier du type .sam"
        
        # Vérifie si le fichier est vide
        if [[ ! -s "$fichier" ]]
         then
            echo "Erreur : Le fichier '$fichier' est vide."
            exit 1
        else
            echo "Le fichier '$fichier' n'est pas vide."
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

            echo "Le fichier '$fichier' respecte le nombre minimum de colonnes."
        fi

            count_reads_in_sam() {
                # Compte le total de reads
                local totalReads=$(awk '!/^@/ {print $2}' "$fichier" | wc -l)
                echo "Le total de reads est: $totalReads"

                # Compte les reads non mappés (flag 4 indique les reads non mappés)
                local unmappedReads=$(awk '$2 == 4 {print $2}' "$fichier" | wc -l)
                echo "Le total de reads non mappés est: $unmappedReads"

                # Compte les reads en duplicité (flag 102 indique les reads en duplicité)
                local duplicatedReads=$(awk '$2 == 102 {print $2}' "$fichier" | wc -l)
                echo "Le total de reads en duplicité est: $duplicatedReads"

                # Compte les reads mappés (les reads non mappés ne sont pas inclus)
                local mappedReads=$(awk '!/^@/ && $2 != 4 {print $2}' "$fichier" | wc -l)
                echo "Le total de reads mappés est: $mappedReads"
            }
            count_reads_in_sam
    else
        echo "Erreur: Ce n'est pas un fichier du type régulier et '.sam'. Insérez un nouveau fichier valide avec les permissions de lecture."
    fi
done

