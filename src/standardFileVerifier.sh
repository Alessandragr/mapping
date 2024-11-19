#!/bin/bash
# -*- coding: utf-8 -*-

# Combien de reads sont mappés par chromosome ? Compter les reads par chromosome. 
# Gestion de plusieurs chromosomes, un dictionnaire pour chaque chromosome
# comment tester

# Où les reads sont-ils mappés ? Est-ce que l’alignement est homogène le long de la séquence de référence ? Compter le nombre de reads par position sur le chromosome. done

# Comment les reads (et paires de reads) sont-ils mappés ? Compter le nombre de reads pour chaque valeur de flag. done

# Avec quelle qualité les reads sont-ils mappés ? Compter les reads pour chaque score de qualité ou par tranche de valeurs de score de mapping. done

# En fonction des questions 2 et 4, implémenter des paramètres pour filtrer les données. Par exemple, ne conserver que les reads complètement mappés ou ceux avec un score de qualité de mapping inférieur à 30.

# Prendre les reads, réaliser le mapping avec un génome de référence et aligner les résultats.

# EXTRA : Projet additionnel : Coder pour enregistrer le **CIGAR** et stocker ces informations afin de connaître le nombre de délétions, insertions, etc., dans chaque séquence.

# EXTRA : Projet additionnel : Analyser la qualité du CIGAR. Plus un CIGAR est long, plus la qualité du read est incertaine, tandis que les CIGAR plus courts indiquent un alignement plus sûr avec le génome de référence.

# EXTRA : Executer Python vs Bash avec Jobs dans un serveur.or

#!/bin/bash

echo "[INFO] Parameters: $*"
for fichier in "$@"; do
    # Check if the filename is valid (alphanumeric, underscores, dashes, and dots allowed)
    if [[ ! "$fichier" =~ ^[a-zA-Z0-9_.-]+$ ]]; then
        echo "[ERROR] Invalid filename: '$fichier'. Ensure it contains only alphanumeric characters, underscores, dashes, or dots."
        continue
    fi

    # Check if the file is in the same directory as the script
    script_dir=$(dirname "$0")
    if [[ ! -f "$script_dir/$fichier" ]]; then
        echo "[ERROR] File '$fichier' not found in the same directory as the script."
        continue
    fi

    # Check if it's a regular file
    if [[ -f "$script_dir/$fichier" ]]; then
        echo "[INFO] The document '$fichier' is a valid file."

        # Check if the file has read permissions
        if [[ -r "$script_dir/$fichier" ]]; then
            echo "[INFO] The file '$fichier' has read permissions."

            # Check if the file has a .sam extension
            if [[ "$fichier" == *.sam ]]; then
                echo "[INFO] The document '$fichier' is of type .sam."

                # Check if the file is empty
                if [[ ! -s "$script_dir/$fichier" ]]; then
                    echo "[ERROR] The file '$fichier' is empty."
                    continue
                else
                    echo "[INFO] The file '$fichier' is not empty."

                    # Validate the number of columns in the first 3 non-header lines
                    line_count=0
                    while IFS= read -r line || [[ -n "$line" ]]; do
                        # Ignore header lines (starting with '@')
                        if [[ $line == @* ]]; then
                            continue
                        fi

                        # Count the number of columns (tab-separated)
                        column_count=$(echo "$line" | awk -F'\t' '{print NF}')

                        # Check if the number of columns is less than 11
                        if (( column_count < 11 )); then
                            echo "[ERROR] Line '$line' in '$fichier' has only $column_count columns. It is not a valid .sam file."
                            continue 2
                        fi

                        # Increment the line counter and stop after the first 3 lines
                        ((line_count++))
                        if (( line_count == 3 )); then
                            break
                        fi
                    done < "$script_dir/$fichier"
                    echo "[SUCCESS] The file '$fichier' meets the minimum column count requirement."
                fi
            else
                echo "[ERROR] The file '$fichier' is not of type '.sam'. Please provide a valid .sam file."
            fi
        else
            echo "[ERROR] The file '$fichier' does not have read permissions. Please provide a file with read permissions."
        fi
    else
        echo "[ERROR] '$fichier' is not a valid file."
    fi
done