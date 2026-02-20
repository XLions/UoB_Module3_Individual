setwd('D:/OneDrive/PROJECTS/UoB/M3/Groupwork');rm(list = ls());gc()
library(tidyverse)
library(readr)
library(progress)

#1.读取数据与预处理数据
if(!file.exists('./Data/expr_mtx.RDS')){
  ##MicroArray
  {
    allProbIntensities<-readr::read_tsv('./Data/DataSet2/allProbIntensities.tsv') %>%
      dplyr::filter(!is.na(GeneSymbols)) %>%
      dplyr::select(-'Reporter.Identifier')
    
    #表达矩阵去重
    array_deduplication<-data.frame(matrix(ncol = ncol(allProbIntensities)-1,
                                           nrow = 0))
    pb <- progress_bar$new(
      total = length(unique(allProbIntensities$GeneSymbols))) # 创建一个进度条对象
    for (i in 1:length(unique(allProbIntensities$GeneSymbols))) {
      array_mtx_temp<-allProbIntensities %>%
        dplyr::filter(GeneSymbols==unique(allProbIntensities$GeneSymbols)[i])
      array_mtx_temp.1<-array_mtx_temp[,-1] %>% colSums()
      array_deduplication<-rbind(
        array_deduplication,
        array_mtx_temp1
      )
      pb$tick()  # 更新进度条
    }
    rownames(array_deduplication)<-unique(allProbIntensities$GeneSymbols)
    colnames(array_deduplication)<-colnames(allProbIntensities)[-1]
    
    rownames(array_deduplication)<-
      str_replace_all(rownames(array_deduplication),'[.]','_')
    rownames(array_deduplication)<-
      str_replace_all(rownames(array_deduplication),'[|]','_')
    rownames(array_deduplication)<-
      str_replace_all(rownames(array_deduplication),'-','_')
  }
  ##RNA-Seq
  {
    log2FPKM<-readr::read_tsv('./Data/DataSet2/log2FPKM.tsv') %>%
      column_to_rownames('00gene_id')
    NewColnames_RNASeq<-c()
    for(i in 1:ncol(log2FPKM)){
      NewColnames_RNASeq<-c(
        NewColnames_RNASeq,
        strsplit(colnames(log2FPKM),'_')[[i]][1]
      )
    }
    colnames(log2FPKM)<-NewColnames_RNASeq
    rownames(log2FPKM)<-
      str_replace_all(rownames(log2FPKM),'[.]','_')
    rownames(log2FPKM)<-
      str_replace_all(rownames(log2FPKM),'[|]','_')
    rownames(log2FPKM)<-
      str_replace_all(rownames(log2FPKM),'-','_')
    rownames(log2FPKM)<-
      str_replace_all(rownames(log2FPKM),'/','_')
  }
  ##保存预处理后表达矩阵
  expr_mtx<-list(
    MicroArray=as.data.frame(array_deduplication),
    RNASeq=as.data.frame(log2FPKM)
  )
  saveRDS(expr_mtx,'./Data/expr_mtx.RDS')
}else{
  expr_mtx<-readRDS('./Data/expr_mtx.RDS')
}
#病人信息
patientInfo_train<-read.delim('./Data/DataSet2/patientInfo_train.tsv')
patientInfo_train$FactorValue..inss.stage.<-
  str_replace_all(
    patientInfo_train$FactorValue..inss.stage.,
    '4S','5'
  ) %>% as.numeric()
patientInfo_test<-readr::read_tsv('./Data/DataSet2/patientInfo_test.tsv')

#2.分析前准备
##2.1.建立因素列名对照表
factor_columnName<-data.frame(
  Factor=c('ID','Sex','Age at Diagnosis','Death from Disease','High Risk',
           'INSS Stage','Progression'),
  ColumnNames=c('ID','FactorValue..Sex.','FactorValue..age.at.diagnosis.',
                'FactorValue..death.from.disease.','FactorValue..high.risk.',
                'FactorValue..inss.stage.','FactorValue..progression.')
)
##2.2.自定义函数
{
  #1).模型评估函数
  modelSummary <- function(real_label, pred_label, prob_pos, positive = NULL) {
    library(tidyverse)
    library(readr)
    library(progress)
    library(caret)
    library(pROC)
    
    real_label <- factor(real_label)
    pred_label <- factor(pred_label, levels = levels(real_label))
    
    if (is.null(positive)) {
      positive <- levels(real_label)[1]
    }
    
    # ---------- AUC ----------
    roc_obj <- pROC::roc(
      response  = real_label,
      predictor = prob_pos,
      levels    = rev(levels(real_label)),
      quiet     = TRUE
    )
    auc_val <- as.numeric(pROC::auc(roc_obj))
    
    # ---------- confusion matrix ----------
    cm <- caret::confusionMatrix(
      data = pred_label,
      reference = real_label,
      positive = positive
    )
    
    acc <- unname(cm$overall["Accuracy"])
    kap <- unname(cm$overall["Kappa"])
    precision <- unname(cm$byClass["Pos Pred Value"])
    recall    <- unname(cm$byClass["Sensitivity"])
    
    f1 <- ifelse(
      (precision + recall) == 0,
      0,
      2 * precision * recall / (precision + recall)
    )
    
    list(
      AUC       = auc_val,
      Accuracy  = acc,
      F1        = f1,
      Precision = precision,
      Recall    = recall,
      Kappa     = kap
    )
  }
  modelSummary_multi <- function(real_label, pred_label, prob_pos) {
    # 加载必要的包
    suppressPackageStartupMessages({
      if (!require(caret)) install.packages("caret"); library(caret)
      if (!require(pROC)) install.packages("pROC"); library(pROC)
      if (!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
      if (!require(tidyr)) install.packages("tidyr"); library(tidyr)
    })
    
    # 确保标签是因子类型，且水平一致
    labels <- factor(c(1,2,3,4,5))
    real_label <- factor(real_label, levels = labels)
    pred_label <- factor(pred_label, levels = labels)
    
    # 1. 混淆矩阵
    confusion_matrix <- confusionMatrix(pred_label, real_label, mode = "everything")
    
    # 2. 整体准确率
    accuracy <- confusion_matrix$overall["Accuracy"]
    
    # 3. 各类别指标
    class_metrics <- confusion_matrix$byClass
    
    # 4. 宏平均和微平均
    macro_avg <- colMeans(class_metrics[, c("Precision", "Recall", "F1")], na.rm = TRUE)
    names(macro_avg) <- paste0("Macro_", names(macro_avg))
    
    # 5. 如果有概率矩阵，计算AUC和ROC
    auc_results <- NULL
    roc_curves <- list()
    
    if (!is.null(prob_pos)) {
      # 确保概率矩阵的列名正确
      if (!all(as.character(labels) %in% colnames(prob_pos))) {
        colnames(prob_pos) <- as.character(labels)
      }
      
      # 计算每个类别的AUC（一对多）
      auc_values <- numeric(length(labels))
      names(auc_values) <- as.character(labels)
      
      for (label in labels) {
        # 创建二分类标签
        binary_real <- as.numeric(real_label == label)
        
        # 获取该类别的预测概率
        prob_class <- prob_pos[, as.character(label)]
        
        # 计算ROC和AUC
        roc_obj <- roc(binary_real, prob_class, quiet = TRUE)
        auc_values[as.character(label)] <- auc(roc_obj)
        roc_curves[[as.character(label)]] <- roc_obj
      }
      
      # 计算宏平均AUC
      macro_auc <- mean(auc_values)
      
      auc_results <- list(
        per_class_auc = auc_values,
        macro_auc = macro_auc,
        roc_curves = roc_curves
      )
    }
    
    # 6. 创建结果汇总
    results <- list(
      confusion_matrix = confusion_matrix$table,
      overall_metrics = c(
        Accuracy = accuracy,
        Kappa = confusion_matrix$overall["Kappa"]
      ),
      class_metrics = class_metrics,
      macro_averages = macro_avg,
      auc_metrics = auc_results
    )
    
    list(
      AUC       = results$auc_metrics$macro_auc
    # Accuracy  = results$overall_metrics["Accuracy"],
    # F1        = results$overall_metrics["Macro_F1"],
    # Precision = results$overall_metrics["Macro_Precision"],
    # Recall    = results$overall_metrics["Macro_Recall"],
    # Kappa     = results$overall_metrics["Kappa"]
    )
  }
  
  #2).选择变量组建训练集表达矩阵和标签的函数
  chooseFactorToBuildTrainsetAndLabel<-
    function(factor,dataset){
      library(tidyverse)
      library(readr)
      library(progress)
      if(dataset=='RNA-Seq'){
        expr<-expr_mtx$RNASeq
      }else if(dataset=='MicroArray'){
        expr<-expr_mtx$MicroArray
      }
      
      library(dplyr)
      
      x<-expr %>% t() %>% as.data.frame() %>%
        rownames_to_column('sample') %>%
        merge(y=patientInfo_train[
          ,c('ID',factor_columnName$ColumnNames[
            which(factor_columnName$Factor==factor)
          ])],by.x='sample',by.y='ID') %>%
        dplyr::rename('FACTOR'=factor_columnName$ColumnNames[
          which(factor_columnName$Factor==factor)
        ])
      y<-factor(x$FACTOR,levels=unique(x$FACTOR))
      x<-x %>% dplyr::select(-'FACTOR') %>%
        column_to_rownames('sample')
      return(
        list(
          x=x,y=y
        )
      )
    }
  #3).使用Lasso和Boruta提取表达矩阵特征
  featureExtract<-function(x,y){
    #01.加载R包与读取数据
    #加载包
    library(tidyverse)
    library(glmnet)
    library(ggsci)
    library(Boruta)
    library(magrittr)
    library(mlbench)
    library(caret)
    library('xgboost')
    library("Matrix")
    library(ImageGP)
    library(ggvenn)
    library(tidyverse)
    library(leaps)
    library(DALEX)
    library(pROC)
    library(ggpubr)
    library(randomForest)
    select=dplyr::select
    
    
    #Models
    {
      #Boruta-------------------------------------------------------------------------
      message('Analyzing - Boruta')
      library(Boruta)
      if(!file.exists('res.Boruta.RDS')){
        message('Analyzing - Buling Boruta models')
        res.Boruta <- Boruta(x=x, y=y, pValue=0.01, mcAdj=T, maxRuns=300)
        message('Analyzing - Saving Boruta models')
        saveRDS(res.Boruta,'res.Boruta.RDS')
      }else{
        res.Boruta <- readRDS('res.Boruta.RDS')
      }
      Boruta <- attStats(res.Boruta) #给出Boruta算法的结果
      boruta_geneids <- Boruta[Boruta$decision=='Confirmed',] %>% rownames(.)
      write_csv(data.frame(gene=boruta_geneids),'01.boruta_res_genelist.csv')
      
    }
    
    message(paste0(
      'Boruta: ',length(boruta_geneids)
    ))
    
    return(list(
      x=x %>% as.data.frame() %>% 
        dplyr::select(boruta_geneids),
      y=y))
  }
  #4).将训练数据分为训练集和验证集
  cutDataToTrainValidation<-function(x,y,seedNumber,valid_percent=0.8){
    set.seed(seedNumber)
    train_sample_n<-ceiling(nrow(x)*valid_percent)
    train_sample_row<-sample(seq(1:nrow(x)),
                             train_sample_n)
    valid_sample_row<-setdiff(seq(1:nrow(x)), train_sample_row)
    return(
      list(
        seed=seedNumber,
        train_mtx=x[train_sample_row,],
        train_label=y[train_sample_row],
        validation_mtx=x[valid_sample_row,],
        validation_label=y[valid_sample_row]
      )
    )
  }
  #5).建立模型的自定义函数
  #5-1).kNN模型函数 - 二分类
  modelBuilding_kNN<-function(train_mtx,
                              train_label,
                              validation_mtx,
                              validation_label){
    library(class)
    library(stringr)
    library(dplyr)
    library(progress)
    library(pROC)
    
    #建立空白的结果评估表格
    Model_result_frame<-data.frame(matrix(nrow = nrow(train_mtx),
                                          ncol = 7))
    colnames(Model_result_frame)<-c('K','AUC','Accuracy','F1','Precision',
                                    'Recall','Kappa')
    
    #建立小循环：测试合适的K值
    for (i in 1:(nrow(train_mtx))) {
      #建立模型
      pred <- knn(train = train_mtx,
                  test = validation_mtx,
                  cl = train_label,
                  k = i,prob = T)
      #结果写入表格
      Model_result_frame[i,]<-
        c(i,
          modelSummary(as.numeric(as.character(validation_label)),
                       as.numeric(as.character(pred)),
                       attr(pred, "prob"),positive = '1'))
    }
    plot_data<-
      pivot_longer(
        Model_result_frame %>%
          dplyr::mutate('Score'=pmin(AUC, Kappa)),
        cols=c(2:8), 
        names_to = "Item", 
        values_to = "Value")
    
    plots_list_ConfusionMtx<-
      ggplot(data = plot_data %>% dplyr::filter(!Item %in% c('AUC','Kappa','Score')))+
      geom_line(mapping = aes(x=K,y=Value))+
      facet_wrap(~Item)+
      scale_y_continuous(limits = c(0,1))+theme_bw()+
      theme(
        axis.title = element_text(face = 'bold')
      )
    plots_list_KappaAUC<-
      ggplot(data = plot_data %>% dplyr::filter(Item %in% c('AUC','Kappa')))+
      geom_line(mapping = aes(x=K,y=Value))+
      facet_wrap(~Item)+
      scale_y_continuous(limits = c(0,1))+theme_bw()+
      theme(
        axis.title = element_text(face = 'bold')
      )
    #确认最佳值
    BestK<-(plot_data %>% dplyr::filter(Item %in% c('Score')))$K[
      which((plot_data %>% dplyr::filter(Item %in% c('Score')))$Value==
              max((plot_data %>% dplyr::filter(Item %in% c('Score')))$Value))
    ]
    if(length(BestK)>1){
      BestK_filter<-Model_result_frame %>% dplyr::filter(K %in% BestK)
      BestK<-BestK_filter$K[which(BestK_filter$AUC==max(BestK_filter$AUC))]
    }
    BestK_AUC<-Model_result_frame$AUC[BestK]
    BestK_Kappa<-Model_result_frame$Kappa[BestK]
    plots_Score<-
      ggplot(data = plot_data %>% dplyr::filter(Item %in% c('Score')))+
      geom_line(mapping = aes(x=K,y=Value))+
      geom_vline(xintercept = BestK,linetype='dashed')+
      annotate("text", x = BestK, y = 1, label = paste0("K = ", BestK)) +
      scale_y_continuous(limits = c(0,1))+theme_bw()+
      labs(y='Socre [min(AUC, Kappa)]')+
      theme(
        axis.title = element_text(face = 'bold')
      )
    result<-
      list(
        train_mtx=train_mtx,
        train_label=train_label,
        BestK=BestK,
        BestQuality=(Model_result_frame %>%
                       dplyr::mutate('Score'=pmin(AUC, Kappa)))[BestK,],
        plots=list(
          CM=plots_list_ConfusionMtx,
          KappaAUC=plots_list_KappaAUC,
          Score=plots_Score
        )
      )
    return(
      result
    )
  }
  #5-2).多项逻辑回归模型函数 - 多分类
  modelBuilding_MultinomialLogistic<-function(train_mtx,
                                              train_label,
                                              validation_mtx,
                                              validation_label,type_data){
    library(glmnet)
    library(nnet)
    library(caret)
    library(pROC)
    
    ## 标准化矩阵
    scale_center <- apply(train_mtx, 2, mean)
    scale_scale  <- apply(train_mtx, 2, sd)
    scale_scale[scale_scale == 0] <- 1
    train_x <- scale(train_mtx, center = scale_center, scale = scale_scale)
    valid_x <- scale(validation_mtx, center = scale_center, scale = scale_scale)
    ## 标签必须是 factor
    train_y <- factor(train_label)
    valid_y <- factor(validation_label)
    
    ## 组合标签与矩阵
    df_train <- data.frame(y = train_y, train_x)
    df_valid <- data.frame(valid_x)
    ## 建立模型
    fit_multinom <- nnet::multinom(
      y ~ .,
      data  = df_train,
      trace = FALSE
    )
    
    ## 预测类别
    pred_class <- predict(fit_multinom, newdata = df_valid)
    ## 预测概率
    pred_prob <- predict(fit_multinom, newdata = df_valid, type = "probs")
    if(type_data=='binary'){
      pred_prob <- as.numeric(pred_prob)
    }
    result<-
      list(
        model=fit_multinom,
        Quality=
          modelSummary(real_label = as.numeric(as.character(valid_y)),
                       pred_label = as.numeric(as.character(pred_class)),
                       prob_pos = pred_prob)
      )
    return(
      result
    )
  }
  #5-3).线性SVM模型函数 - 多分类
  modelBuilding_LinearSVM<-function(train_mtx,
                                    train_label,
                                    validation_mtx,
                                    validation_label,type_data){
    library(e1071)
    library(caret)
    library(pROC)
    
    
    ## 标准化矩阵
    scale_center <- apply(train_mtx, 2, mean)
    scale_scale  <- apply(train_mtx, 2, sd)
    scale_scale[scale_scale == 0] <- 1
    train_x <- scale(train_mtx, center = scale_center, scale = scale_scale)
    valid_x <- scale(validation_mtx, center = scale_center, scale = scale_scale)
    ## 标签必须是 factor
    train_y <- factor(train_label)
    valid_y <- factor(validation_label)
    
    ## 建立模型
    svm_fit <- svm(
      x = train_x,
      y = train_y,
      kernel = "linear",
      cost = 1,              # 可调：0.1 / 1 / 10
      scale = FALSE,
      probability = TRUE     # 为 AUC 做准备
    )
    
    ## 预测类别
    pred_class <- predict(svm_fit, newdata = valid_x)
    ## 预测概率
    pred_obj <- predict(
      svm_fit,
      newdata = valid_x,
      probability = TRUE
    )
    pred_prob <- attr(pred_obj, "probabilities")
    pred_prob <- pred_prob[, levels(train_y)] #确保列顺序与因子水平一致
    if(type_data=='binary'){
      pred_prob <- as.numeric(pred_prob[,1])
    }
    result<-
      list(
        model=svm_fit,
        features=colnames(train_mtx),
        Quality=
          modelSummary(real_label = as.numeric(as.character(valid_y)),
                       pred_label = as.numeric(as.character(pred_class)),
                       prob_pos = pred_prob)
      )
    return(
      result
    )
  }
  #5-4).随机森林 - 多分类
  modelBuilding_RF<-function(train_mtx,
                             train_label,
                             validation_mtx,
                             validation_label,
                             seedNumber,type_data){
    library(e1071)
    library(caret)
    library(pROC)
    
    set.seed(seedNumber)
    
    ## 标准化矩阵
    scale_center <- apply(train_mtx, 2, mean)
    scale_scale  <- apply(train_mtx, 2, sd)
    scale_scale[scale_scale == 0] <- 1
    train_x <- scale(train_mtx, center = scale_center, scale = scale_scale)
    valid_x <- scale(validation_mtx, center = scale_center, scale = scale_scale)
    ## 标签必须是 factor
    train_y <- factor(train_label)
    valid_y <- factor(validation_label)
    
    library(randomForest)
    
    set.seed(123)
    
    rf_fit <- randomForest(
      x = train_x,
      y = train_y,
      ntree = 1000,
      importance = TRUE
    )
    #预测分类和概率
    pred_class <- predict(rf_fit, newdata = valid_x)
    pred_prob  <- predict(rf_fit, newdata = valid_x, type = "prob")
    if(type_data=='binary'){
      pred_prob <- as.numeric(pred_prob[,1])
    }
    #特征重要性
    imp <- importance(rf_fit)
    imp_df <- data.frame(
      Gene = rownames(imp),
      MeanDecreaseGini = imp[, "MeanDecreaseGini"]
    )
    imp_df <- imp_df[order(-imp_df$MeanDecreaseGini), ]
    top_n <- 20
    imp_top <- head(imp_df, top_n)
    imp_plot<-
      ggplot(imp_top,
             aes(x = reorder(Gene, MeanDecreaseGini),
                 y = MeanDecreaseGini)) +
      geom_col(fill = "#E64B35") +
      coord_flip() +
      labs(
        title = "Random Forest Feature Importance",
        x = NULL,
        y = "Mean Decrease in Gini"
      ) +
      theme_bw()+theme(
        plot.title = element_text(face = 'bold', hjust=0.5),
        axis.title.x = element_text(face = 'bold')
      )
    
    result<-
      list(
        model=rf_fit,
        features=colnames(train_mtx),
        Quality=
          modelSummary(real_label = as.numeric(as.character(valid_y)),
                       pred_label = as.numeric(as.character(pred_class)),
                       prob_pos = pred_prob),
        imp=list(
          table=importance(rf_fit),
          plot=imp_plot
        )
      )
    return(
      result
    )
  } 
  
  #6).摘取AUC最大模型的预测结果
  takePredByMaxAUCModel<-function(path){
    path0<-getwd()
    setwd(path)
    data1<-read.csv('ModelPerformancs.csv') %>% 
      dplyr::filter(Item=='AUC')
    model<-data1$Model[which(data1$Value==max(data1$Value))]
    model1<-ifelse(model=='kNN','kNN',
                   ifelse(model=='ML','MultinomialLogistic',
                          ifelse(model=='LS','LinearSVM','RF')))
    data2<-read.csv('PredResult.csv')
    setwd(path0)
    return(
      data2 %>% 
        dplyr::select(c('ID',model1)) %>%
        dplyr::rename('Pred'=model1)
    )
  }
}
#2.3.组建测试集表达矩阵
testExpr<-list(
  RNASeq=expr_mtx$RNASeq %>%
    dplyr::select(patientInfo_test$ID) %>% t() %>% as.data.frame(),
  MicroArray=expr_mtx$MicroArray %>%
    dplyr::select(patientInfo_test$ID) %>% t() %>% as.data.frame()
)
saveRDS(testExpr,'testExpr.RDS')


#3.1.分析：主循环-特征提取
for (i1 in c('RNA-Seq','MicroArray')) {
  message(paste0('Prepare Data - ',i1,' All'))
  if(!dir.exists(i1)){dir.create(i1)};setwd(i1)
  #准备完整测试集
  if(i1=='RNA-Seq'){
    test_Expr<-testExpr$RNASeq
  }else if(i1=='MicroArray'){
    test_Expr<-testExpr$MicroArray
  }
  for (i2 in 4:nrow(factor_columnName)) {
    message(paste0('Prepare Data - ',i1,' ',
                   c('Death','HighRisk','Stage','Progression')[i2-3]))
    #设置目录
    dir_names<-c('Death','HighRisk','Stage','Progression')[i2-3]
    if(!dir.exists(dir_names)){dir.create(dir_names)};setwd(dir_names)
    #根据时间戳设置种子
    if(!file.exists('seedNumber_Time.txt')){
      seedNumber_Time<-(as.numeric(Sys.time())+i2)
      write(seedNumber_Time,'seedNumber_Time.txt')
    }else{
      seedNumber_Time<-as.numeric(readLines('./seedNumber_Time.txt'))
    }
    #生成数据
    DataRaw<-
      chooseFactorToBuildTrainsetAndLabel(
        factor = factor_columnName$Factor[i2],
        dataset = i1
      )
    #提取特征
    if(!dir.exists('FearureExtract')){dir.create('FearureExtract')};setwd('FearureExtract')
    Data0<-featureExtract(DataRaw$x,DataRaw$y)
    saveRDS(Data0,'Data0.RDS')
    setwd('../')
    setwd('../')
  }
  setwd('../')
}
#特征数绘图
{
  ggplot(data=data.frame(X=rep(c('Death','HighRisk','Stage','Progression'),2),
                         Y=c(65,137,70,51,59,106,54,34),
                         Z=rep(c('RNA-Seq','MicroArray'),each=4)),
         mapping = aes(x=X,y=Y,fill = Z))+
    geom_bar(stat = 'identity', 
             #柱状图位置并排:
             position = position_dodge(width=0.9), #组内柱子间隔
             width = 0.7,      #设置柱子宽度,使变量之间分开
             color='black')+ 
    scale_fill_manual(values = c('indianred','lightskyblue'))+
    theme_bw()+
    labs(x='Factor',y='Number of features',fill='Dataset')+
    theme(axis.title = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold', hjust = 0.5))
  ggsave('FeatureExtracting.pdf',height = 5,width = 6)
  ggsave('FeatureExtracting.png',height = 5,width = 6, dpi=600)
}

#3.2.分析 - 二分类因素
{
  #3.2.1.分析：主循环-二分类因素的模型建立
  for (i1 in c('RNA-Seq','MicroArray')) {
    message(paste0('Prepare Data - ',i1,' All'))
    if(!dir.exists(i1)){dir.create(i1)};setwd(i1)
    #准备完整测试集
    if(i1=='RNA-Seq'){
      test_Expr<-testExpr$RNASeq
    }else if(i1=='MicroArray'){
      test_Expr<-testExpr$MicroArray
    }
    for (i2 in 4:(nrow(factor_columnName)-1)) {
      message(paste0('Prepare Data - ',i1,' ',
                     c('Death','HighRisk','Progression')[i2-3]))
      #设置目录
      dir_names<-c('Death','HighRisk','Progression')[i2-3]
      if(!dir.exists(dir_names)){dir.create(dir_names)};setwd(dir_names)
      #根据时间戳设置种子
      if(!file.exists('seedNumber_Time.txt')){
        seedNumber_Time<-(as.numeric(Sys.time())+i2)
        write(seedNumber_Time,'seedNumber_Time.txt')
      }else{
        seedNumber_Time<-as.numeric(readLines('./seedNumber_Time.txt'))
      }
      #生成数据
      DataRaw<-
        chooseFactorToBuildTrainsetAndLabel(
          factor = factor_columnName$Factor[i2],
          dataset = i1
        )
      #读取提取特征结果
      Data0<-readRDS('./FearureExtract/Data0.RDS')
      FeatureGeneTrain<-colnames(Data0$x)
      #随机分割为训练集和验证集
      Data1<-
        cutDataToTrainValidation(
          Data0$x,Data0$y,
          seedNumber=seedNumber_Time,
          valid_percent=0.8)
      saveRDS(Data1,'TrainValidation.RDS') #保存
      
      #训练模型
      message(paste0('Model kNN - ',i1,' ',
                     c('Death','HighRisk','Progression')[i2-3]))
      Model_kNN<-
        modelBuilding_kNN(train_mtx = Data1$train_mtx,
                          train_label = Data1$train_label,
                          validation_mtx = Data1$validation_mtx,
                          validation_label = Data1$validation_label)
      saveRDS(Model_kNN,'Model_kNN.RDS') 
      
      message(paste0('Model MLogic - ',i1,' ',
                     c('Death','HighRisk','Progression')[i2-3]))
      Model_MultinomialLogistic<-
        modelBuilding_MultinomialLogistic(train_mtx = Data1$train_mtx,
                                          train_label = Data1$train_label,
                                          validation_mtx = Data1$validation_mtx,
                                          validation_label = Data1$validation_label,
                                          type_data='binary')
      saveRDS(Model_MultinomialLogistic,'Model_MultinomialLogistic.RDS')
      
      message(paste0('Model SVM - ',i1,' ',
                     c('Death','HighRisk','Progression')[i2-3]))
      Model_LinearSVM<-
        modelBuilding_LinearSVM(train_mtx = Data1$train_mtx,
                                train_label = Data1$train_label,
                                validation_mtx = Data1$validation_mtx,
                                validation_label = Data1$validation_label,
                                type_data='binary')
      saveRDS(Model_LinearSVM,'Model_LinearSVM.RDS')
      
      message(paste0('Model RF - ',i1,' ',
                     c('Death','HighRisk','Progression')[i2-3]))
      Model_RF<-
        modelBuilding_RF(train_mtx = Data1$train_mtx,
                         train_label = Data1$train_label,
                         validation_mtx = Data1$validation_mtx,
                         validation_label = Data1$validation_label,
                         seedNumber = seedNumber_Time,
                         type_data='binary')
      saveRDS(Model_RF,'Model_RF.RDS')
      
      #模型绘图
      ModelsPlot_Data<-data.frame(
        Model=rep(c('kNN','ML','LS','RF'),each=6),
        Item=rep(c('AUC','Accuracy','F1','Precision',
                   'Recall','Kappa'),4),
        Value=c(
          as.numeric(distinct(Model_kNN$BestQuality[2:7])),
          as.numeric(Model_MultinomialLogistic$Quality),
          as.numeric(Model_LinearSVM$Quality),
          as.numeric(Model_RF$Quality)
        )
      )
      write.csv(ModelsPlot_Data,'ModelPerformancs.csv',row.names = F)
      performances_plot<-
        ggplot(data=ModelsPlot_Data,
               mapping = aes(x=Item,y=Value,fill=Model))+
        geom_bar(stat = 'identity', 
                 #柱状图位置并排:
                 position = position_dodge(width=0.9), #组内柱子间隔
                 width = 0.7,      #设置柱子宽度,使变量之间分开
                 color='black')+ 
        scale_fill_manual(values=c('gold','lightskyblue','#3CB371','#CD5C5C'))+
        geom_text(aes(label=round(Value,3)),size=4,
                  position = position_dodge(width = 0.9), #相应的注释宽度也调整
                  vjust=-0.3)+    #调节注释高度    
        theme_bw()+  
        theme(axis.title = element_text(colour = 'black'),
              legend.title = element_text(face = 'bold',hjust=0.5))
      ggsave('performances_plot.pdf',height=6,width=13,performances_plot)
      
      setwd('../')
    }
    setwd('../')
  }
  #3.2.2.分析：主循环-二分类因素的预测
  for (i1 in c('RNA-Seq','MicroArray')) {
    message(paste0('Prepare Data - ',i1,' All'))
    if(!dir.exists(i1)){dir.create(i1)};setwd(i1)
    #准备完整测试集
    if(i1=='RNA-Seq'){
      test_Expr<-testExpr$RNASeq
    }else if(i1=='MicroArray'){
      test_Expr<-testExpr$MicroArray
    }
    for (i2 in 4:(nrow(factor_columnName)-1)) {
      message(paste0('Prepare Data - ',i1,' ',
                     c('Death','HighRisk','Progression')[i2-3]))
      #设置目录
      dir_names<-c('Death','HighRisk','Progression')[i2-3]
      if(!dir.exists(dir_names)){dir.create(dir_names)};setwd(dir_names)
      #根据时间戳设置种子
      if(!file.exists('seedNumber_Time.txt')){
        seedNumber_Time<-(as.numeric(Sys.time())+i2)
        write(seedNumber_Time,'seedNumber_Time.txt')
      }else{
        seedNumber_Time<-as.numeric(readLines('./seedNumber_Time.txt'))
      }
      
      #读取模型
      message('Reading Models')
      Model_kNN<-readRDS('Model_kNN.RDS')
      Model_MultinomialLogistic<-readRDS('Model_MultinomialLogistic.RDS')
      Model_LinearSVM<-readRDS('Model_LinearSVM.RDS')
      Model_RF<-readRDS('Model_RF.RDS')
      
      #预测
      message(paste0('Predicting Testsets - ',i1,' ',
                     c('Death','HighRisk','Progression')[i2-3]))
      PredResult<-data.frame(
        ID=patientInfo_test$ID,
        kNN=NA,MultinomialLogistic=NA,LinearSVM=NA,RF=NA
      )
      test_Expr1<-test_Expr
      test_Expr2<-test_Expr
      PredResult$kNN<-
        knn(train = Model_kNN$train_mtx,
            test = test_Expr1 %>% dplyr::select(colnames(Model_kNN$train_mtx)),
            cl = Model_kNN$train_label,
            k = Model_kNN$BestK[1],prob = T) %>% as.character() %>% as.numeric() 
      PredResult$MultinomialLogistic<-
        predict(Model_MultinomialLogistic$model, 
                newdata = test_Expr2 %>% dplyr::select(Model_MultinomialLogistic[["model"]][["coefnames"]][-1]))
      PredResult$LinearSVM<-
        predict(Model_LinearSVM$model, 
                newdata = test_Expr1 %>% dplyr::select(Model_LinearSVM$features))
      PredResult$RF<-
        predict(Model_RF$model, 
                newdata = test_Expr1 %>% dplyr::select(Model_RF$features))
      saveRDS(PredResult,'PredResult.RDS')
      write.csv(PredResult,'PredResult.csv',row.names = F)
      
      save.image('Image.RData')
      setwd('../')
    }
    setwd('../')
  }
}

#3.3.分析 - 多分类因素
{
  #3.2.1.分析：主循环-多分类因素的模型建立
  for (i1 in c('RNA-Seq','MicroArray')) {
    message(paste0('Prepare Data - ',i1,' All'))
    if(!dir.exists(i1)){dir.create(i1)};setwd(i1)
    #准备完整测试集
    if(i1=='RNA-Seq'){
      test_Expr<-testExpr$RNASeq
    }else if(i1=='MicroArray'){
      test_Expr<-testExpr$MicroArray
    }
    message(paste0('Prepare Data - ',i1,' ',
                   'Stage'))
    #设置目录
    dir_names<-'Stage'
    if(!dir.exists(dir_names)){dir.create(dir_names)};setwd(dir_names)
    #根据时间戳设置种子
    if(!file.exists('seedNumber_Time.txt')){
      seedNumber_Time<-(as.numeric(Sys.time())+i2)
      write(seedNumber_Time,'seedNumber_Time.txt')
    }else{
      seedNumber_Time<-as.numeric(readLines('./seedNumber_Time.txt'))
    }
    #生成数据
    DataRaw<-
      chooseFactorToBuildTrainsetAndLabel(
        factor = factor_columnName$Factor[i2],
        dataset = i1
      )
    #读取提取特征结果
    Data0<-readRDS('./FearureExtract/Data0.RDS')
    FeatureGeneTrain<-colnames(Data0$x)
    #随机分割为训练集和验证集
    Data1<-
      cutDataToTrainValidation(
        Data0$x,Data0$y,
        seedNumber=seedNumber_Time,
        valid_percent=0.8)
    saveRDS(Data1,'TrainValidation.RDS') #保存
    
    #训练模型
    train_mtx = Data1$train_mtx
    train_label = Data1$train_label
    validation_mtx = Data1$validation_mtx
    validation_label = Data1$validation_label
    
    message(paste0('Model MLogic - ',i1,' ',
                   'Stage'))
    #Model_MultinomialLogistic
    {
      library(glmnet)
      library(nnet)
      library(caret)
      library(pROC)
      
      ## 标准化矩阵
      scale_center <- apply(train_mtx, 2, mean)
      scale_scale  <- apply(train_mtx, 2, sd)
      scale_scale[scale_scale == 0] <- 1
      train_x <- scale(train_mtx, center = scale_center, scale = scale_scale)
      valid_x <- scale(validation_mtx, center = scale_center, scale = scale_scale)
      ## 标签必须是 factor
      train_y <- factor(train_label)
      valid_y <- factor(validation_label)
      
      ## 组合标签与矩阵
      df_train <- data.frame(y = train_y, train_x)
      df_valid <- data.frame(valid_x)
      ## 建立模型
      fit_multinom <- nnet::multinom(
        y ~ .,
        data  = df_train,
        trace = FALSE
      )
      
      ## 预测类别
      pred_class <- predict(fit_multinom, newdata = df_valid)
      ## 预测概率
      pred_prob <- predict(fit_multinom, newdata = df_valid, type = "probs")
      
      Model_MultinomialLogistic<-
        list(
          model=fit_multinom,
          Quality=
            modelSummary_multi(real_label = as.numeric(as.character(valid_y)),
                               pred_label = as.numeric(as.character(pred_class)),
                               prob_pos = pred_prob)
        )
      saveRDS(Model_MultinomialLogistic,'Model_MultinomialLogistic.RDS')
    }
    
    message(paste0('Model SVM - ',i1,' ',
                   'Stage'))
    #SVM
    {
      library(e1071)
      library(caret)
      library(pROC)
      
      
      ## 标准化矩阵
      scale_center <- apply(train_mtx, 2, mean)
      scale_scale  <- apply(train_mtx, 2, sd)
      scale_scale[scale_scale == 0] <- 1
      train_x <- scale(train_mtx, center = scale_center, scale = scale_scale)
      valid_x <- scale(validation_mtx, center = scale_center, scale = scale_scale)
      ## 标签必须是 factor
      train_y <- factor(train_label)
      valid_y <- factor(validation_label)
      
      ## 建立模型
      svm_fit <- svm(
        x = train_x,
        y = train_y,
        kernel = "linear",
        cost = 1,              # 可调：0.1 / 1 / 10
        scale = FALSE,
        probability = TRUE     # 为 AUC 做准备
      )
      
      ## 预测类别
      pred_class <- predict(svm_fit, newdata = valid_x)
      ## 预测概率
      pred_obj <- predict(
        svm_fit,
        newdata = valid_x,
        probability = TRUE
      )
      pred_prob <- attr(pred_obj, "probabilities")
      pred_prob <- pred_prob[, levels(train_y)] #确保列顺序与因子水平一致
      
      Model_LinearSVM<-
        list(
          model=svm_fit,
          features=colnames(train_mtx),
          Quality=
            modelSummary_multi(real_label = as.numeric(as.character(valid_y)),
                               pred_label = as.numeric(as.character(pred_class)),
                               prob_pos = pred_prob)
        )
      saveRDS(Model_LinearSVM,'Model_LinearSVM.RDS')
    }
    
    message(paste0('Model RF - ',i1,' ',
                   'Stage'))
    {
      library(e1071)
      library(caret)
      library(pROC)
      
      
      ## 标准化矩阵
      scale_center <- apply(train_mtx, 2, mean)
      scale_scale  <- apply(train_mtx, 2, sd)
      scale_scale[scale_scale == 0] <- 1
      train_x <- scale(train_mtx, center = scale_center, scale = scale_scale)
      valid_x <- scale(validation_mtx, center = scale_center, scale = scale_scale)
      ## 标签必须是 factor
      train_y <- factor(train_label)
      valid_y <- factor(validation_label)
      
      library(randomForest)
      
      set.seed(123)
      
      rf_fit <- randomForest(
        x = train_x,
        y = train_y,
        ntree = 1000,
        importance = TRUE
      )
      #预测分类和概率
      pred_class <- predict(rf_fit, newdata = valid_x)
      pred_prob  <- predict(rf_fit, newdata = valid_x, type = "prob")
      #特征重要性
      imp <- importance(rf_fit)
      imp_df <- data.frame(
        Gene = rownames(imp),
        MeanDecreaseGini = imp[, "MeanDecreaseGini"]
      )
      imp_df <- imp_df[order(-imp_df$MeanDecreaseGini), ]
      top_n <- 20
      imp_top <- head(imp_df, top_n)
      imp_plot<-
        ggplot(imp_top,
               aes(x = reorder(Gene, MeanDecreaseGini),
                   y = MeanDecreaseGini)) +
        geom_col(fill = "#E64B35") +
        coord_flip() +
        labs(
          title = "Random Forest Feature Importance",
          x = NULL,
          y = "Mean Decrease in Gini"
        ) +
        theme_bw()+theme(
          plot.title = element_text(face = 'bold', hjust=0.5),
          axis.title.x = element_text(face = 'bold')
        )
      
      Model_RF<-
        list(
          model=rf_fit,
          features=colnames(train_mtx),
          Quality=
            modelSummary_multi(real_label = as.numeric(as.character(valid_y)),
                               pred_label = as.numeric(as.character(pred_class)),
                               prob_pos = pred_prob),
          imp=list(
            table=importance(rf_fit),
            plot=imp_plot
          )
        )
      saveRDS(Model_RF,'Model_RF.RDS')
    }
    
    #模型绘图
    ModelsPlot_Data<-data.frame(
      Model=rep(c('ML','LS','RF'),1),
      Item=rep(c('AUC'),3),
      Value=c(
        as.numeric(Model_MultinomialLogistic$Quality),
        as.numeric(Model_LinearSVM$Quality),
        as.numeric(Model_RF$Quality)
      )
    )
    write.csv(ModelsPlot_Data,'ModelPerformancs.csv',row.names = F)
    performances_plot<-
      ggplot(data=ModelsPlot_Data,
             mapping = aes(x=Item,y=Value,fill=Model))+
      geom_bar(stat = 'identity', 
               #柱状图位置并排:
               position = position_dodge(width=0.9), #组内柱子间隔
               width = 0.7,      #设置柱子宽度,使变量之间分开
               color='black')+ 
      scale_fill_manual(values=c('gold','lightskyblue','#3CB371','#CD5C5C'))+
      geom_text(aes(label=round(Value,3)),size=4,
                position = position_dodge(width = 0.9), #相应的注释宽度也调整
                vjust=-0.3)+    #调节注释高度    
      theme_bw()+  
      theme(axis.title = element_text(colour = 'black'),
            legend.title = element_text(face = 'bold',hjust=0.5))
    ggsave('performances_plot.pdf',height=6,width=13,performances_plot)
    
    setwd('../')
    setwd('../')
  }
  #3.2.2.分析：主循环-多分类因素的预测
  for (i1 in c('RNA-Seq','MicroArray')) {
    message(paste0('Prepare Data - ',i1,' All'))
    if(!dir.exists(i1)){dir.create(i1)};setwd(i1)
    #准备完整测试集
    if(i1=='RNA-Seq'){
      test_Expr<-testExpr$RNASeq
    }else if(i1=='MicroArray'){
      test_Expr<-testExpr$MicroArray
    }
    message(paste0('Prepare Data - ',i1,' ',
                   'Stage'))
    #设置目录
    dir_names<-'Stage'
    if(!dir.exists(dir_names)){dir.create(dir_names)};setwd(dir_names)
    #根据时间戳设置种子
    if(!file.exists('seedNumber_Time.txt')){
      seedNumber_Time<-(as.numeric(Sys.time())+i2)
      write(seedNumber_Time,'seedNumber_Time.txt')
    }else{
      seedNumber_Time<-as.numeric(readLines('./seedNumber_Time.txt'))
    }
    
    #读取模型
    message('Reading Models')
    Model_MultinomialLogistic<-readRDS('Model_MultinomialLogistic.RDS')
    Model_LinearSVM<-readRDS('Model_LinearSVM.RDS')
    Model_RF<-readRDS('Model_RF.RDS')
    
    #预测
    message(paste0('Predicting Testsets - ',i1,' ',
                   'Stage'))
    PredResult<-data.frame(
      ID=patientInfo_test$ID,
      kNN=NA,MultinomialLogistic=NA,LinearSVM=NA,RF=NA
    )
    test_Expr1<-test_Expr
    test_Expr2<-test_Expr
    PredResult$MultinomialLogistic<-
      predict(Model_MultinomialLogistic$model, 
              newdata = test_Expr2 %>% dplyr::select(Model_MultinomialLogistic[["model"]][["coefnames"]][-1]))
    PredResult$LinearSVM<-
      predict(Model_LinearSVM$model, 
              newdata = test_Expr1 %>% dplyr::select(Model_LinearSVM$features))
    PredResult$RF<-
      predict(Model_RF$model, 
              newdata = test_Expr1 %>% dplyr::select(Model_RF$features))
    saveRDS(PredResult,'PredResult.RDS')
    write.csv(PredResult,'PredResult.csv',row.names = F)
    
    save.image('Image.RData')
    setwd('../')
    setwd('../')
  }
}

#4.统计结果：选择AUC最大的模型
for(i1 in c('RNA-Seq','MicroArray')){
  pred_df<-data.frame(ID=patientInfo_test$ID)
  for (i2 in c('Death','HighRisk','Progression','Stage')) {
    pred_i2<-takePredByMaxAUCModel(path = paste0('./',i1,'/',i2))
    colnames(pred_i2)[2]<-'Pred'
    pred_df<-merge(x=pred_df,y=pred_i2,by='ID')
    colnames(pred_df)[ncol(pred_df)]<-i2
  }
  write.csv(pred_df,paste0('Pred_',i1,'.csv'))
}

