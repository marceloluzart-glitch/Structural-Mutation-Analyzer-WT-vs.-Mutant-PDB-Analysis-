# Structural Mutation Analyzer (WT vs. Mutant PDB Analysis)
# Analisador Estrutural de Mutações (Análise PDB WT vs. Mutante)

> **Confidential Research Repository / Repositório de Pesquisa Confidencial:** Automated computational pipeline for structural comparison, global/local RMSD calculation, residue-wise displacement profiling, and 3D visualization.  
> *Pipeline computacional automatizado para comparação estrutural, cálculo de RMSD global/local, perfil de deslocamento por resíduo e visualização 3D.*

---

## 🚀 Overview / Visão Geral

* **[EN]** When analyzing protein variants (such as point mutations modeled via AlphaFold), quantifying structural deviations accurately is crucial to understanding functional impacts. This pipeline performs automated comparative analysis between a **Wild-Type (WT)** protein structure and its **Mutant counterpart**.  
*Note: Raw structural coordinates and unpublished proprietary datasets are strictly excluded from this repository to protect ongoing research.*

* **[PT]** Ao analisar variantes proteicas (como mutações pontuais modeladas via AlphaFold), quantificar os desvios estruturais com precisão é fundamental para compreender impactos funcionais. Este pipeline realiza análises comparativas automatizadas entre a estrutura de uma proteína **Selvagem (WT)** e sua **versão mutante**.  
*Nota: Coordenadas estruturais brutas e conjuntos de dados proprietários não publicados são estritamente excluídos deste repositório para resguardar a pesquisa em andamento.*

### Key Features / Principais Funcionalidades:
* **Global & Local Superposition / Superposição Global e Local:** Computes global RMSD across all aligned $C_\alpha$ atoms and local RMSD within a customizable sliding window around the mutation hotspot. *(Calcula o RMSD global em todos os átomos $C_\alpha$ alinhados e o RMSD local em uma janela deslizante ao redor do sítio da mutação).*
* **Residue-wise Distance Profiling / Perfil de Distância por Resíduo:** Calculates positional Euclidean distances for every amino acid residue to map structural shifts. *(Calcula distâncias euclidianas posicionais para cada aminoácido para mapear alterações conformacionais).*
* **Automated Plotting / Geração de Gráficos:** Generates distance profile plots highlighting structural perturbations near the mutation site. *(Produz perfis de distância destacando perturbações estruturais próximas ao ponto de mutação).*
* **Interactive 3D Visualization / Visualização 3D Interativa:** Renders interactive molecular structures for structural inspection using `py3Dmol`. *(Renderiza estruturas moleculares interativas utilizando `py3Dmol`).*

---

## 🛠️ Tech Stack / Tecnologias Utilizadas

* **Python 3.x**
* **BioPython** (PDB parsing and structural alignment via `Superimposer` / *Parsing de PDB e alinhamento estrutural*)
* **NumPy** (Vectorized coordinate handling and Euclidean distance math / *Manipulação vetorizada de coordenadas*)
* **Matplotlib** (Structural displacement profiling / *Mapeamento de deslocamento estrutural*)
* **py3Dmol** (Interactive 3D molecular rendering / *Renderização molecular 3D interativa*)

---

## 📋 Installation & Requirements / Instalação e Requisitos

Clone the repository and install the dependencies:  
*Clone o repositório e instale as dependências:*

```bash
git clone [https://github.com/seu-usuario/structural-mutation-analyzer.git](https://github.com/seu-usuario/structural-mutation-analyzer.git)
cd structural-mutation-analyzer
pip install biopython numpy matplotlib py3Dmol
