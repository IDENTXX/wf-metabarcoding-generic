wf-metabarcoding-generic-v2
Ein Nextflow-basierter Workflow für EPI2ME Desktop zur Analyse von Oomyceten (COX-Gen, ITS) und Fusarium (TEF1). Der Workflow wurde speziell für die Nutzung unter Windows mit WSL2 und Docker optimiert.

Features
Hybrider Modus: Dateioperationen laufen lokal (WSL), Bioinformatik-Tools laufen in Docker-Containern, um Pfadprobleme unter Windows zu vermeiden.

Intelligente Taxonomie: Automatische Bereinigung von GenBank-Accession-Nummern und korrekte Darstellung von Varietäten sowie forma specialis (z. B. f. sp. lycopersici).

Erweiterbar: Einfaches Hinzufügen neuer Datenbanken über den Katalog in der main.nf.

Voraussetzungen
Docker Desktop: Version 4.62.0 oder neuer (mit aktivierter WSL2-Integration).

EPI2ME Desktop: Zur grafischen Ausführung des Workflows.

Speicherort: Datenbanken sollten unter D:/Epi2Me_Datenbanken/ liegen.

Installation & Setup
Klone dieses Repository in deinen lokalen EPI2ME-Workflow-Ordner.

Stelle sicher, dass deine Referenzdatenbanken im FASTA-Format vorliegen:

oomycetes_cox_ref.fasta

fusarium_TEF_ref.fasta

oomycetes_its_ref.fasta

Nutzung
Öffne EPI2ME Desktop.

Wähle diesen Workflow aus.

Wähle die gewünschte Datenbank im Dropdown-Menü aus.

Gib den Pfad zu deinen Nanopore-Reads (.fastq.gz) an.

Klicke auf Launch.

Output
Die Ergebnisse befinden sich im Ordner output/summary:

wf-metagenomics-counts-species.csv: Die Haupttabelle mit bereinigten Artnamen und Read-Counts pro Barcode.

abundance_matrix.csv: Eine kompakte Matrix für statistische Auswertungen.
