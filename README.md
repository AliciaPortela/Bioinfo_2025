# Parte 1: Tests de selección

## 1. Introducción

Esta práctica comienza a integrar datos genómicos en su contexto
ambiental mediante la inferencia de la huella de la selección natural en
el genoma de distintas poblaciones de humanos. Estos tests de selección
funcionan descomponiendo la distancia genética entre poblaciones en dos
partes que llamaremos alfa y beta.

**Beta** es el componente poblacional, es la fracción de distancia
genética que se atribuye a la mera existencia de las poblaciones sobre
la que se mide. La divergencia o convergencia genética entre poblaciones
se puede deber a muchos factores que nada tienen que ver con la
selección natural. Este componente beta se infiere buscando el modelo
demográfico (cada test de selección tiene en cuenta distintos
parámetros) que mejor explica las distancias genéticas observadas entre
poblaciones. **Alfa** es el componente de locus y se traduce en magnitud
del efecto de la selección. Lo que hacen los tests de selección es
utilizar beta como hipótesis nula, calcular distancia genética para cada
locus en el genoma y atribuir desvíos de esa beta a eventos de
selección. Algunos tests incluso discriminan entre desvíos positivos y
negativos (más o menos distancia genética que la esperada en función de
beta, respectivamente), de manera que infieren selección divergente o
purificadora.

### ¿Qué tests vamos a usar?

Vamos a aprender a utilizar dos aproximaciones muy distintas para
inferir selección a lo largo del genoma: controlando el efecto de la
genealogía ([GRoSS](https://github.com/FerRacimo/GRoSS)) y el de la
estructura poblacional jerárquica de nuestra muestra
([Flink](https://bitbucket.org/wegmannlab/flink/wiki/browse/)).

**GRoSS** infiere el componente beta del que hablábamos antes midiendo
distancias entre poblaciones a través de su genealogía. Lo que hace es
calcular unas frecuencias alélicas esperadas en cada población en
función de la deriva genética que se va acumulando a lo largo de su
historia con una chi-cuadrado. Atribuye los desvíos de esas frecuencias
a la selección. Además, infiere en qué rama de la genealogía es más
probable que se haya producido el evento de selección mediante la
reconstrucción de estados ancestrales de dicho desvío.

**Flink** infiere el componente beta montando un modelo demográfico
valiéndose del ligamiento que presentan los distintos sitios de genoma
de la muestra. Flink no utiliza todos los sitios del genoma como si
fueran independientes porque no lo son, muchos de esos sitios segregan
juntos y ofrecen información perfectamente redundante. Los modelos
demográficos se obtienen mediante inferencia bayesiana a partir de las
distribuciones de múltiples parámetros importantes, como tasas de
migración entre las distintas poblaciones, tasas de mutación,
desequilibrio en el tamaño de las poblaciones, etc.

## 2. Preparación del entorno

Arranca el ordenador y elige la partición de disco de Ubuntu, vamos a
trabajar casi todo en bash hoy.

Abre la terminal y muévete con *cd* hasta la ubicación de trabajo que
elijas que vas a utilizar.

    # Clonar el repositorio del curso
    git clone https://github.com/AliciaPortela/Bioinfo_2025

    # Entrar al directorio de trabajo
    cd Bioinfo_2025

Familiarízate con la estructura del repositorio y descarga el software
necesario para correr los tests de selección.

    # Muévete a la carpeta desde donde vas a correr GRoSS y descárgalo
    cd 2.GRoSS
    git clone https://github.com/FerRacimo/GRoSS

    # Haz lo mismo con flink
    cd ../../3.Flink
    git clone --recurse-submodules https://bitbucket.org/WegmannLab/flink.git

    # y compílalo (conviértelo en ejecutable)
    cd flink
    make clean
    make

## 3. Ejecución de GRoSS

GRoSS se ejecuta a través de un script en R (`GRoSS.R`) que espera como
entrada tres argumentos:

-   `-e`: archivo con los conteos de distintos genotipos por población.
-   `-r`: archivo con el árbol de relaciones entre poblaciones en
    formato `.graph`.
-   `-o`: nombre del archivo de salida.

Los archivos necesarios ya están disponibles en la carpeta
`2.GRoSS/inputs/`. Puedes encontrar información sobre cómo conseguir el
input de este software
[aquí](https://github.com/FerRacimo/GRoSS/blob/master/VCFtoGRoSS.md)

    # Ejecutar GRoSS (ajustando nombres y rutas de archivos)
    Rscript GRoSS.R -e [ruta a los conteos] -r [ruta al graph] -o [ruta al output].tsv

Nota que especificamos una extensión concreta para el output y la razón
es que es la que se esperan los scripts de ploteo de los resultados.

## 4. Ejecución de Flink

Flink tiene distintas tareas disponibles entre las que está la de
estimar efecto de la selección sobre el grado de divergencia o
convergencia para cada uno de los sitios del genoma (`estimate`). Puedes
encontrar información sobre cómo conseguir el formato de input del
programa
[aquí](https://bitbucket.org/wegmannlab/flink/wiki/Prepare%20Input%20Data).

    $dir_flink/flink \
      task=estimate \
      numThreads=10 \
      B="(-2.0,1.8)" \
      A_max=4.0 \
      lnMu_g="(-4.0,0.0)" \
      lnNu_g="(-5.0,0.0)" \
      s_max=14 \
      beta="(-2.0,1.8)" \
      alpha_max=4.0 \
      numIterations=500000 \
      burnin=300000 \
      thinning=100 \
      lnkappa="(-10.0,-0.1)" \
      logFile=logfile_human \
      sigmaProp_mu=0.005 \
      sigmaProp_nu=0.05 \
      sigmaProp_kappa=0.05 \
      data=input.flink.txt

[Aquí](https://bitbucket.org/wegmannlab/flink/wiki/Parameter%20Estimation)
puedes consultar el significado de cada uno de los parámetros.

## 5. Resultados

Hay muchas cosas que podemos hacer con los resultados de estos tests
como comprobar si algún sitio se corresponde con alguna mutación no
homóloga dentro de un gen con función reconocida
([BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html),
[Uniprot](https://www.uniprot.org/)), o correr análisis de asociación
con fenotipos para intentar interpretar esas funciones
([GWAS](https://www.cog-genomics.org/plink/2.0/tutorials/qwas1)).

Lo que vamos a hacer nosotros en la práctica va a ser comparar la
estructura genética que obtenemos con variación neutral con la que
obtenemos utilizando sitios del genoma que inferimos que están sujetos a
selección. De esta manera podemos interpretar distintas fuentes de
divergencia o convergencia entre las poblaciones, incluidos procesos de
adaptación local.

### Estructura genética basada en variación neutral

En la carpeta `1.Genomic_data` del repositorio del curso podéis
encontrar la base de datos original con toda la diversidad genética del
trozo de genoma que estamos utilizando: `human_data_QC.vcf`.

Hemos corrido un PCA para que podáis visualizarlo directamente en R,
pero os dejamos aquí una breve explicación de cómo hacerlo para que la
tengáis a mano:

    # Descarga y deja listo plink para correr
    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250615.zip
    unzip plink_linux_x86_64_20250615.zip

    # Corre el PCA
    cd plink_linux_x86_64_20250615
    ./plink --vcf [ruta al vcf] --double-id --set-missing-var-ids @:# \
    --make-bed --pca --out [prefijo del output]

El argumento `--make-bed` es para conseguir un formato de archivo que no
es necesario para conseguir el PCA, pero resulta bastante útil a la hora
de hacer otros análisis de estructura que no vamos a ver en este módulo,
como [Admixture](https://dalexander.github.io/admixture/). El output de
este análisis está en la misma carpeta y se llama `PCA_data.txt`. Este
archivo, junto con las correspondencias de cada individuo con sus
poblaciones (`popfile.txt`), es lo que necesitáis para visualizar el
PCA.

    library(ggplot2)
    library(dplyr)

    setwd("")
    # Leer tabla de PCA
    pca_data <- read.table("PCA_data.txt", header=TRUE, stringsAsFactors = FALSE)

    # Leer tabla de poblaciones
    pop_data <- read.table("popfile.txt", sep = ",", header=TRUE, stringsAsFactors = FALSE)

    # Unir tablas por la columna (IID)
    pca_pop <- left_join(pca_data, pop_data, by = "IID")

    # Plot PC1 vs PC2 por continente
    plot <- ggplot(pca_pop, aes(x=PC1, y=PC2, color=CONTINENT)) +
      geom_point(size=3) +
      theme_minimal() +
      labs(title="PCA_all", x="PC1", y="PC2") +
      theme(legend.title = element_text(size=12),
            legend.text = element_text(size=10)) + stat_ellipse()

<img src="https://raw.githubusercontent.com/AliciaPortela/Bioinfo_2025/main/1.Genomic_data/PCA_all.png" width="60%" />

### Estructura genética basada en SNPs sujetos a selección (GRoSS)

Para extraer los SNPs que están sujetos a selección según GRoSS, hemos
preparado un script en R que podéis encontrar en
`2.GRoSS/Outputs/extract_outliers.R`. Lo que hace este script es
básicamente encontrar el mínimo valor de significación para el efecto de
la selección en las frecuencias alélicas de cada SNP (recuerda que se
realiza el análisis una vez por rama de la genealogía). Después, extrae
una lista de SNPs cuya significación mínima sea menor a un valor que
pretende ser lo menos arbitrario posible y que se estima a partir de un
screeplot con todos los valores de significación. Una vez obtenida esa
lista de SNPs (`SNPs_under_sel.txt`), la base de datos original es
fácilmente filtrable utilizando
[bcftools](https://samtools.github.io/bcftools/howtos/index.html) (no
hace falta que lo hagáis porque ya está todo listo para visualizar el
PCA, que también está ya hecho, de los SNPs sujetos a selección según
GRoSS):

    # comprime y crea el índice del vcf original
    bgzip -c [ruta al vcf] > [ruta al vcf].gz
    tabix -p vcf [ruta al vcf].gz

    ./bcftools view -R SNPs_under_sel.txt -o ./human_data_SNPs_under_sel.vcf -O v ./hgdp_tgp_sgdp_chr12_p.dated_287.vcf.gz

Una vez con la base de datos filtrada, simplemente habría que usar plink
para obtener el PCA y R para plotearlo como hemos hecho antes.

<img src="https://raw.githubusercontent.com/AliciaPortela/Bioinfo_2025/main/2.GRoSS/Outputs/PCA_SNPs_under_sel.png" width="60%" />

Además, en la carpeta `2.GRoSS` tenéis un script para plotear la
significación del efecto de la selección sobre las frecuencias alélicas
de cada SNP y para cada rama de la genealogía: `plot_outliers.R`. Lo que
hace este script es básicamente plotear cada uno de los p\_values del
efecto de la selección sobre las frecuencias alélicas de cada SNP para
cada una de las ramas de la genealogía.

    library("ggplot2")
    library("reshape")

    setwd("")

    table <- read.table("human.tsv",header=TRUE)
    table <- as.data.frame(table)

    finaltab <- table

    CHR <- finaltab[,1]
    START <- finaltab[,2]
    END <- finaltab[,3]
    MIDPOINT <- (finaltab[,2]+finaltab[,3])/2
    newtab <- finaltab[,seq(4,dim(finaltab)[2])]
    newtab <- newtab[,seq(dim(newtab)[2]/2+1,dim(newtab)[2])]
    POS <- seq(1,dim(newtab)[1])
    melttab <- melt(newtab)
    melttab[,2] <- -log10(melttab[,2])
    melttab <- cbind(melttab, POS,CHR,START,END)
    melttab[,1] <- factor(melttab[,1])
    names(melttab) <- c("Branches","Pvals","SNPID","CHR","START","END")
    totalbranches <- length(unique(melttab[,1]))

    plot1 <- ggplot(melttab, aes(x=SNPID, y=Pvals, fill=Branches, shape=Branches,colour=Branches))+
    geom_point()+
    scale_shape_manual(values=rep(c(4,15,17,18,19),100)[seq(1,totalbranches)])
    print(plot1)
    ggsave("outliers.png", plot = plot1, width = 12, height = 10, dpi = 300)

Para interpretarlo correctamente, tenéis que plotear también el .graph
que habíais usado como entrada del análisis y así ver a qué ramas
corresponden los valores de significación más bajos:

    library(igraph)
    library(ggplot2)

    setwd("") 

    lines <- readLines("human.graph")

    # Extrae las aristas (líneas que comienzan con 'edge')
    edges <- grep("^edge", lines, value = TRUE)
    edge_df <- do.call(rbind, strsplit(edges, "\t"))
    edge_df <- edge_df[, 3:4]  # Solo columnas de origen y destino

    # Crea el grafo dirigido
    g <- graph_from_edgelist(as.matrix(edge_df), directed = TRUE)

    plot(g,
         layout = layout_as_tree(g),
         vertex.label.cex = 0.8,
         vertex.size = 15,
         vertex.color = "lightblue",
         edge.arrow.size = 0.4)

<img src="https://raw.githubusercontent.com/AliciaPortela/Bioinfo_2025/main/2.GRoSS/Outputs/outliers.png" height="300px">

<img src="https://raw.githubusercontent.com/AliciaPortela/Bioinfo_2025/main/2.GRoSS/Outputs/graph_plot.png" height="300px">

### Estructura genética basada en SNPs sujetos a selección (Flink)

En este caso, los autores publican un script para extraer los SNPs que
Flink infiere que están bajo selección, lo tenéis en
`3.Flink/Outputs/findSel.R`. Lo puedes correr en bash y lo único que
requiere son los siguientes parámetros:

    chr="chr12" 
    group="" # Hierarchy para el nivel más alto, revisa el input para el resto de grupos y nómbralo aquí
    index=0 # Numera el grupo de antes según el orden del input, de 0-6
    sel="" # Tipo de selección: divergent or balancing
    prob=0.99999 # la probabilidad de descubrimiento falso (1-significación)

    Rscript findSel.R $chr $group $index $sel $prob

De esta manera se obtiene una lista de SNPs con la que filtrar la base
de datos original y poder correr un PCA como ya hemos visto.

<img src="https://raw.githubusercontent.com/AliciaPortela/Bioinfo_2025/main/3.Flink/Outputs/PCA_flink_balanced.png" width="45%">
<img src="https://raw.githubusercontent.com/AliciaPortela/Bioinfo_2025/main/3.Flink/Outputs/PCA_flink_divergent.png" width="45%">

Los autores también publican un script para visualizar la intensidad de
cada tipo de selección por grupos, lo tenéis en
`3.Flink/Outputs/plot_loci_under_sel.R`. Lo hemos modificado aquí
ligeramente para que añada alguna información relevante:

    setwd("")

    args <- commandArgs(trailingOnly = TRUE);
    name <- "12";
    rm(args);

    library(data.table)

    freq_pos<-read.table("inputFlink.txt",header=F,skip = 2)
    data_flink_g<-fread("Posterior_alphas_group_0.txt",header=T) # Cambia esto para plotear distintos alphas

    data_split <- split(freq_pos, freq_pos$V1)

    pdf(paste("plotSel",name,".pdf",sep=""),width = 14, height = 3.6);
    layout(matrix(c(1,1,1,2,2,2,3,3,3), 9,1, byrow = TRUE))
    par(mar=c(1,6.5,1.0,2.0),mgp=c(3.5,1.5,0),oma=c(6,0,0,0))

    numloci<-nrow(as.data.frame(data_split[name]))
    divPostProb<-vector(mode="numeric",length=numloci)
    balPostProb<-vector(mode="numeric",length=numloci)

    alphastart_d <- read.table("start_alphas.txt", skip = 1, header = F)
    alphastart <- alphastart_d[,2]
    names(alphastart) = alphastart_d[,1]

    calcDivSel <- function(data,numloci,divPostProb) {
      rows<-nrow(data)
      divPostProb=(colSums(data[,alphastart[name]:(alphastart[name]+numloci)]>0))/(rows+1)
    }

    calcBalSel <- function(data,numloci,balPostProb) {
      rows<-nrow(data)
      balPostProb=(colSums(data[,alphastart[name]:(alphastart[name]+numloci)]<0))/(rows+1)
    }

    plotdata <- function(freq,divPostProb,balPostProb,namepop) {
      plot(-log(1.0-divPostProb,base=10),main=,ylab="",xlab="",cex.main=2.5,las=1, type='l',ylim=c(0,3), xlim=c(0,5000),cex.lab=2.5,cex.axis=2,xaxs="i",yaxs="i",col="orange2",xaxt="n",yaxt="n")
      lines(-log(1.0-balPostProb,base=10),col="dodgerblue")
      legend("topright",legend=namepop, bty ="n", pch=NA,cex=1.6,inset=c(-0.0001,-0.15)) 
      mtext(expression(paste(-log[10]~q)), side=2, line=4,cex=1.3)
      axis(2, seq(0,3,2),las=2,cex.axis=2)
      abline(h=2.0,lty="dashed") 
    }

    divPostProb<-calcDivSel(data_flink_g,numloci,divPostProb)
    balPostProb<-calcBalSel(data_flink_g,numloci,balPostProb)
    plotdata(freq_pos[name]$V2,divPostProb,balPostProb,name)
    axis(1, las=TRUE, cex.axis=2.0)
    mtext(paste("Locus",sep=""), side=1, line=4,cex=1.6,outer=T)

    plotdata <- function(freq,divPostProb,balPostProb,namepop) {
      plot(-log(1.0-divPostProb,base=10),main=,ylab="",xlab="",cex.main=2.5,las=1, type='l',ylim=c(0,3), xlim=c(5000,10000),cex.lab=2.5,cex.axis=2,xaxs="i",yaxs="i",col="orange2",xaxt="n",yaxt="n")
      lines(-log(1.0-balPostProb,base=10),col="dodgerblue")
      legend("topright",legend=namepop, bty ="n", pch=NA,cex=1.6,inset=c(-0.0001,-0.15)) 
      mtext(expression(paste(-log[10]~q)), side=2, line=4,cex=1.3)
      axis(2, seq(0,3,2),las=2,cex.axis=2)
      abline(h=2.0,lty="dashed") 
    }

    plotdata(freq_pos[name]$V2,divPostProb,balPostProb,name)
    axis(1, las=TRUE, cex.axis=2.0)
    mtext(paste("Locus",sep=""), side=1, line=4,cex=1.6,outer=T)

    dev.off();

[Ver PDF del resultado de
selección](https://github.com/AliciaPortela/Bioinfo_2025/raw/main/3.Flink/Outputs/plotSelNA.pdf)
