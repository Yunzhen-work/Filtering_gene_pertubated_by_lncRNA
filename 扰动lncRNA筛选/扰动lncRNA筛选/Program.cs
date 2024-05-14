using System;
using System.Collections;
using System.Collections.Generic;
using System.Data;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using RDotNet;

//2019.08.11 计算完了rewire score和100次随机，还有p值，但bonferroni-corrected的p值R返回的都是整数

namespace 扰动lncRNA筛选
{
    class Program
    {
        static string workingDirectory = @"E:\BRCA_workspace_0\";

        //如果有第三步保存的文件就填写这个，没有的话就改成空，也就是连续两个引号
        static string step3OutputFilePath = "";
        //static string step3OutputFilePath = workingDirectory + "step3_filtered_triplet.txt";

        static string lncRNAExpressionFilePath = workingDirectory + "BRCA-lncRNA-purity.txt";
        static string microRNAExpressionFilePath = workingDirectory + "BRCA-microRNA-purity.txt";
        static string mRNAExpressionFilePath = workingDirectory + "BRCA-mRNA-purity.txt";
        static string tripletFilePath = workingDirectory + "triplet.txt";
        static string outputPath3 = workingDirectory + "step3_filtered_triplet.txt";
        static string outputPath4 = workingDirectory + "step4_filtered_triplet.txt";
        

        static char[] splitter1 = new char[] {'\t'}, splitter2 = new char[] {'\t', '-'};
        
        static DataTable RNAFilesReader(string filePath, string tableName, Hashtable requiredList = null,
            bool assignPrimaryKey = true)
        {
            StreamReader reader = null;
            DataTable dt = new DataTable(tableName);
            string row_read = null;
            string[] row_split = null;
            bool notFinishReading = true;

            while (notFinishReading)
            {
                try
                {
                    reader = new StreamReader(filePath);
                    Console.WriteLine("正在读取" + tableName + "..." + filePath);
                    if (requiredList != null)
                    {
                        Console.WriteLine("只读取triplet中有的项...");
                    }

                    //header
                    row_read = reader.ReadLine();
                    row_split = row_read.Split(splitter1, StringSplitOptions.RemoveEmptyEntries);
                    foreach (string columnName in row_split)
                    {
                        dt.Columns.Add(columnName);
                    }

                    //context
                    while (!reader.EndOfStream)
                    {
                        row_read = reader.ReadLine();
                        row_split = row_read.Split(splitter1, StringSplitOptions.RemoveEmptyEntries);
                        if (requiredList != null)
                        {
                            if (requiredList.Contains(row_split[0]))
                                dt.Rows.Add(row_split);
                        }
                        else
                        {
                            dt.Rows.Add(row_split);
                        }
                    }

                    //finished
                    reader.Close();
                    notFinishReading = false;
                    Console.WriteLine("读取完毕，共" + dt.Rows.Count + "行（不含表头），" + dt.Columns.Count + "列\n");

                    if (assignPrimaryKey)
                        dt.PrimaryKey = new DataColumn[] {dt.Columns[0]};
                }
                catch (Exception e)
                {
                    reader?.Close();
                    Console.WriteLine(e.Message);
                    Console.WriteLine("文件读取错误，请输入" + tableName + "文件的正确路径：");
                    filePath = Console.ReadLine();
                    notFinishReading = true;
                    //throw;
                }

            }

            return dt;
        }

        static DataTable RNAFilterByIQR(DataTable dt, double IQRThreshold)
        {
            DataTable filteredDt = dt.Clone();
            Console.Write("正在筛选" + dt.TableName + "，条件IQR > " + IQRThreshold + " ...");
            foreach (DataRow dataRow in dt.Rows)
            {
                List<double> expressionArray = new List<double>();
                for (int i = 1; i < dataRow.ItemArray.Length; i++)
                {
                    expressionArray.Add(Convert.ToDouble(dataRow.ItemArray[i]));
                }

                expressionArray.Sort();
                int sampleAmount = expressionArray.Count;
                double IQR = expressionArray[sampleAmount - 1] - expressionArray[Convert.ToInt32(sampleAmount * 0.25)];
                if (IQR > IQRThreshold)
                {
                    filteredDt.Rows.Add(dataRow.ItemArray);
                }
            }

            Console.WriteLine("筛选完后剩余" + filteredDt.Rows.Count + "行");
            return filteredDt;
        }

        static Hashtable ConvertDataTableToHashtable(DataTable dt, int colNo)
        {
            Hashtable ht = new Hashtable();
            foreach (DataRow dataRow in dt.Rows)
            {
                if (!ht.Contains(dataRow[colNo].ToString()))
                    ht.Add(dataRow[colNo].ToString(), 1);
            }

            return ht;
        }

        static Hashtable ConvertDataTableToHashtable(DataTable dt, string colName)
        {
            Hashtable ht = new Hashtable();
            foreach (DataRow dataRow in dt.Rows)
            {
                if (!ht.Contains(dataRow[colName].ToString()))
                    ht.Add(dataRow[colName].ToString(), 1);
            }

            return ht;
        }

        static string IsDifferentiallyExpressedByLimmaInR(REngine rEngine, string nameOfTestedRNA,
            DataTable RNAexpressionTable, string highGroupSampleString, string lowGroupSampleString)
        {
            DataRow[] testedRNAExpressionArray =
                RNAexpressionTable.Select(RNAexpressionTable.TableName + " = '" + nameOfTestedRNA + "'");
            DataRow testedRNAExpression = null;
            if (testedRNAExpressionArray.Length == 0)
                return "FILTERED";
            else
                testedRNAExpression = testedRNAExpressionArray[0];

            //System.IO.Path.GetTempPath()
            StreamWriter divideWriter = null, sampleWriter = null;

            divideWriter =
                new StreamWriter(System.IO.Path.GetTempPath() +
                                 "\\divide.txt"); //输出divide.txt，两【列】：sample，category，先high后low
            sampleWriter =
                new StreamWriter(System.IO.Path.GetTempPath() +
                                 "\\sample.txt"); //输出sample.txt，两【行】：sample，(testedRNA Name)，先high后low

            divideWriter.WriteLine("sample\tcategory");
            string rowHeader = "sample";
            string rowValues = nameOfTestedRNA;


            string[] row_split = highGroupSampleString.Split(splitter1, StringSplitOptions.RemoveEmptyEntries);
            foreach (string sampleName in row_split)
            {
                divideWriter.WriteLine(sampleName + "\thigh");
                rowHeader += "\t" + sampleName;
                rowValues += "\t" + testedRNAExpression[sampleName];
            }

            row_split = lowGroupSampleString.Split(splitter1, StringSplitOptions.RemoveEmptyEntries);
            foreach (string sampleName in row_split)
            {
                divideWriter.WriteLine(sampleName + "\tlow");
                rowHeader += "\t" + sampleName;
                rowValues += "\t" + testedRNAExpression[sampleName];
            }

            divideWriter.Close();
            sampleWriter.WriteLine(rowHeader);
            sampleWriter.WriteLine(rowValues);
            sampleWriter.Close();

            rEngine.Evaluate("exprSet = read.table(\"" + System.IO.Path.GetTempPath().Replace(@"\", @"\\") +
                             "sample.txt\", sep = \"\\t\", header = TRUE)").AsDataFrame();
            rEngine.Evaluate("rownames(exprSet)=exprSet[,1]");
            rEngine.Evaluate("exprSet=exprSet[,-1]");
            rEngine.Evaluate("library(limma)");
            rEngine.Evaluate("data1 = read.table(\"" + System.IO.Path.GetTempPath().Replace(@"\", @"\\") +
                             "divide.txt\", sep = \"\\t\", header = TRUE)").AsDataFrame();
            rEngine.Evaluate("design=model.matrix(~factor(data1$category))");
            rEngine.Evaluate("fit=lmFit(exprSet,design)");
            rEngine.Evaluate("fit=eBayes(fit)");
            rEngine.Evaluate("options(digits=4)");
            var output = rEngine.Evaluate("topTable(fit,coef=2,adjust='BH')").AsDataFrame();

            //若返回结果logFC < -0.58表明是差异的
            double logFC = (double) output[0][0];
            if ((logFC < -0.58))
                return "YES";
            else
                return "NO";
        }

        /// <summary>
        /// 计算SCC值，[0]为SCC值，[1]为pValue
        /// </summary>
        /// <param name="rEngine"></param>
        /// <param name="nameOfMicroRNA"></param>
        /// <param name="nameOfMRNA"></param>
        /// <param name="microRNAexpressionTable"></param>
        /// <param name="mRNAexpressionTable"></param>
        /// <param name="sampleString"></param>
        /// <returns></returns>
        static double[] CalculateSCCValueByR(REngine rEngine, string nameOfMicroRNA, string nameOfMRNA,
            DataTable microRNAexpressionTable, DataTable mRNAexpressionTable, string sampleString)
        {

            DataRow microRNAExpressionRow =
                microRNAexpressionTable.Select(microRNAexpressionTable.TableName + " = '" + nameOfMicroRNA + "'")[0];

            DataRow mRNAExpressionRow =
                mRNAexpressionTable.Select(mRNAexpressionTable.TableName + " = '" + nameOfMRNA + "'")[0];


            string[] sampleNames = sampleString.Split(splitter1, StringSplitOptions.RemoveEmptyEntries);
            double[] microRNAExpressionArray = new double[sampleNames.Length];
            double[] mRNAExpressionArray = new double[sampleNames.Length];

            for (int i = 0; i < sampleNames.Length; i++)
            {
                microRNAExpressionArray[i] = Convert.ToDouble(microRNAExpressionRow[sampleNames[i]]);
                mRNAExpressionArray[i] = Convert.ToDouble(mRNAExpressionRow[sampleNames[i]]);
            }

            NumericVector x = new NumericVector(rEngine, microRNAExpressionArray);
            NumericVector y = new NumericVector(rEngine, mRNAExpressionArray);
            rEngine.SetSymbol("x", x);
            rEngine.SetSymbol("y", y);
            GenericVector result = rEngine.Evaluate("cor.test(x, y, method=\"spearman\")").AsList();
            double pValue = result[2].AsNumeric()[0]; //pValue是第3个值，下标[2]
            double SCC = result[3].AsNumeric()[0]; //相关系数是第4个值，下标[3]

            return new double[] {SCC, pValue};
        }

        static double FisherTransformation(double Rvalue)
        {
            return 0.5 * Math.Log(((1 + Rvalue) / (1 - Rvalue)), Math.E);
        }

        static double PValueFromNormalGaussianDistribution(REngine rEngine, double Miu, bool isXSmallerThanMiu, bool isTwoTailedTest)
        {
            NumericVector x = new NumericVector(rEngine, new double[] {Miu});
            rEngine.SetSymbol("x", x);
            NumericVector y = rEngine.Evaluate("pnorm(x, mean = 0, sd = 1)").AsNumeric();
            if (isXSmallerThanMiu)
            {
                return isTwoTailedTest ? -1 : y[0];
            }
            else
            {
                return isTwoTailedTest ? (1 - y[0]) * 2 : 1 - y[0];
            }
        }

        static void Main(string[] args)
        {
            StreamWriter writer = null;
            //之前做了microRNA-mRNA在整体样本中的相关性筛选，接下来做identification of lncRNA modulators across cancers

            //string step3OutputFilePath = "";
            //！！所有文件都要有表头！lncRNA、mRNA、microRNA都要有，大小写要一致

            #region 读取数据


            DataTable triplet = RNAFilesReader(tripletFilePath, "triplet", null, false);
            Hashtable lncRNAinTriplet = ConvertDataTableToHashtable(triplet, "lncRNA");
            Hashtable microRNAinTriplet = ConvertDataTableToHashtable(triplet, "microRNA");
            Hashtable mRNAinTriplet = ConvertDataTableToHashtable(triplet, "mRNA");

            //只读取triplet中有的
            DataTable lncRNA = RNAFilesReader(lncRNAExpressionFilePath, "lncRNA", lncRNAinTriplet);
            DataTable microRNA = RNAFilesReader(microRNAExpressionFilePath, "microRNA", microRNAinTriplet);
            DataTable mRNA = RNAFilesReader(mRNAExpressionFilePath, "mRNA", mRNAinTriplet);

            #endregion

            //第一步：筛选lncRNA、microRNA和mRNA，log2(IQR) > 0.58 => IQR > 1.495

            #region

            /*
            DataTable filtered_lncRNA = lncRNA;
            DataTable filtered_microRNA = microRNA;
            DataTable filtered_mRNA = mRNA;
            */

            DataTable filtered_lncRNA = RNAFilterByIQR(lncRNA, 1.495);
            DataTable filtered_microRNA = RNAFilterByIQR(microRNA, 1.495);
            DataTable filtered_mRNA = RNAFilterByIQR(mRNA, 1.495);

            #endregion

            //第二步：根据每个lncRNA的表达值排序样本，得出前25%和后25%的样本列表顺序（H-group，L-group），改用Hashtable存放！间隔符为\t
            //目前只存放SampleName，如果之后需要Expression值的话每个group就再追加一个Hashtable来存放

            #region

            Console.Write("根据每个lncRNA的表达值排序样本，得出前25%和后25%的样本列表顺序（H-group，L-group）...");

            int oneForthQuantile = Convert.ToInt32(lncRNA.Columns.Count * 0.25); //也是每个group的元素数量

            Hashtable highGroup = new Hashtable(); //存的都是SampleName，只有25%的量，用hashtable存放吧
            Hashtable lowGroup = new Hashtable();
            
            foreach (DataRow dataRow in filtered_lncRNA.Rows)
            {
                //行列转置，然后按表达值数据列降序排序
                DataTable dtSort = new DataTable();
                dtSort.Columns.Add("Sample");
                dtSort.Columns.Add("Expression", typeof(double));
                for (int i = 1;
                    i < lncRNA.Columns.Count;
                    i++) //lncRNA和filtered的两个datatable列结构是一样的，可以直接互用；第0列是lncRNA的Name，不放入表内
                {
                    dtSort.Rows.Add(lncRNA.Columns[i].ColumnName, Convert.ToDouble(dataRow.ItemArray[i]));
                }

                dtSort.DefaultView.Sort = "Expression DESC"; //多列排序的话条件用逗号分隔
                dtSort = dtSort.DefaultView.ToTable();

                //存储highGroup和lowGroup到Hashtable中
                int lGroupOffset = (lncRNA.Columns.Count - 1) - oneForthQuantile;
                string highGroupSampleString = dtSort.Rows[0]["Sample"].ToString();//接sampleString，先放上第0行
                string lowGroupSampleString = dtSort.Rows[lGroupOffset]["Sample"].ToString();

                for (int i = 1; i < oneForthQuantile; i++)//接着上面的第0行
                {
                    highGroupSampleString += "\t" + dtSort.Rows[i]["Sample"];
                    lowGroupSampleString += "\t" + dtSort.Rows[i + lGroupOffset]["Sample"]; //下标位置靠后，所以平移一下
                }

                highGroup.Add(dataRow.ItemArray[0].ToString(), highGroupSampleString);
                lowGroup.Add(dataRow.ItemArray[0].ToString(), lowGroupSampleString);

                //调试，输出排序前后的Sample Name
                /*
                for (int i = 1; i < lncRNA.Columns.Count; i++)
                {
                    Console.WriteLine(lncRNA.Columns[i].ColumnName + "    " + dataRow.ItemArray[i] + "   " +
                                      dtSort.Rows[i - 1][0] + " " + dtSort.Rows[i - 1][1]);
                }
                Console.WriteLine();
                */

            }

            Console.WriteLine("完毕");
            Console.WriteLine();

            #endregion

            //第三步：从H-group L-group列表和triplet列表筛选，要求microRNA在两组内没有差异表达，但mRNA在两组内要有差异表达（P-adjusted < 0.01 and fold change > 1.5）
            //之后用R里的一个包，先按对应的h-group和l-group里的Sample顺序输出出来
            Console.WriteLine("从H-group L-group列表和triplet列表筛选，要求microRNA在两组内没有差异表达，但mRNA在两组内要有差异表达...");

            //初始化RDotNet
            REngine.SetEnvironmentVariables();
            REngine rEngine = REngine.GetInstance();
            rEngine.Initialize();
            rEngine.Evaluate("setwd(\"" +
                             System.IO.Path.GetTempPath().TrimEnd(new char[] { '\\' }).Replace(@"\", @"\\") + "\")");


            DataTable filteredByStep3_triplet = triplet.Clone();
            int n = 0;

            if (step3OutputFilePath == "")
            {
                //microRNA无差异表达，mRNA有差异表达。遇到不符合要求的triplet就不放到新datatable里
                //有无差异表达要调用R里的limma包

                
                foreach (DataRow tripletRow in triplet.Rows)
                {
                    n++;
                    Console.Write("第" + n + "个triplet     ");
                    Console.SetCursorPosition(0, Console.CursorTop);
                    string current_microRNA_name = tripletRow["microRNA"].ToString();
                    string current_mRNA_name = tripletRow["mRNA"].ToString();
                    string current_lncRNA_name = tripletRow["lncRNA"].ToString();

                    if (highGroup[current_lncRNA_name] == null)
                        continue;

                    string IsMicroRNADifferentiallyExpressed = IsDifferentiallyExpressedByLimmaInR(rEngine,
                        current_microRNA_name, filtered_microRNA,
                        highGroup[current_lncRNA_name].ToString(), lowGroup[current_lncRNA_name].ToString());

                    //microRNA要求无差异表达，若有差异则不要，或已被IQR筛选掉了也不要
                    if ((IsMicroRNADifferentiallyExpressed == "YES") ||
                        (IsMicroRNADifferentiallyExpressed == "FILTERED"))
                        continue;

                    string IsMRNADifferentiallyExpressed = IsDifferentiallyExpressedByLimmaInR(rEngine,
                        current_mRNA_name, filtered_mRNA,
                        highGroup[current_lncRNA_name].ToString(), lowGroup[current_lncRNA_name].ToString());

                    //mRNA要求有差异表达，若无差异则不要，或已被IQR筛选掉了也不要
                    if ((IsMRNADifferentiallyExpressed == "NO") || (IsMRNADifferentiallyExpressed == "FILTERED"))
                        continue;

                    //满足条件，放入新的datatable里
                    filteredByStep3_triplet.Rows.Add(tripletRow.ItemArray);

                }

                Console.WriteLine("完毕！第3步筛选后的triplet共" + filteredByStep3_triplet.Rows.Count + "行");
                Console.WriteLine();

                writer = new StreamWriter(outputPath3);
                writer.WriteLine(triplet.Columns[0].ColumnName + "\t" + triplet.Columns[1].ColumnName + "\t" +
                                 triplet.Columns[2].ColumnName);
                foreach (DataRow dataRow in filteredByStep3_triplet.Rows)
                {
                    writer.WriteLine(string.Join("\t", dataRow.ItemArray));
                }

                writer.Close();
            }

            else
            {
                Console.WriteLine("直接读取上次输出的文件..."+ step3OutputFilePath);
                filteredByStep3_triplet = RNAFilesReader(step3OutputFilePath, "triplet", null, false);
            }
            //第四步：对筛选出的triplet测试是否lncRNA有扰动。分别计算H-group L-group里对应的每一对microRNA-mRNA的Spearman系数Rhigh和Rlow，
            //要求任意一个是强相关（|R| > 某个threshold，默认0.2）并且两个必须都是负相关（因已有研究成果表明)(R < 0)，
            //并且lncRNA表达高的时候（high组）二者相关性变弱，也就是|Rlow|-|Rhigh|>某个threshold，默认0.3
            Console.WriteLine("计算R_high，R_low，rewire score并随机排列100次计算rewire score...");

            writer = new StreamWriter(outputPath4);

            //Header
            writer.Write(triplet.Columns[0].ColumnName + "\t" + triplet.Columns[1].ColumnName + "\t" +
                             triplet.Columns[2].ColumnName + "\tR_high\tP(R_high)\tR_low\tP(R_low)\tP Condition\tSCC Condition");
            //Rewire值也要输出
            writer.WriteLine("\tFhigh\tFlow\tRewire score\tP(Rewire score)\tP-adjusted(Bonferroni-corrected)");


            n = 0;
            //DataTable filteredByStep4_triplet = triplet.Clone();
            //filteredByStep4_triplet.Columns.Add("P(Rewire score)");
            List<double> pValues = new List<double>();
            foreach (DataRow tripletRow in filteredByStep3_triplet.Rows)
            {
                n++;
                Console.Write("第" + n + "个triplet     ");
                Console.SetCursorPosition(0, Console.CursorTop);
                string current_microRNA_name = tripletRow["microRNA"].ToString();
                string current_mRNA_name = tripletRow["mRNA"].ToString();
                string current_lncRNA_name = tripletRow["lncRNA"].ToString();

                writer.Write(string.Join("\t", tripletRow.ItemArray));

                double[] Rhigh = CalculateSCCValueByR(rEngine, current_microRNA_name, current_mRNA_name,
                    filtered_microRNA, filtered_mRNA, highGroup[current_lncRNA_name].ToString());

                double[] Rlow = CalculateSCCValueByR(rEngine, current_microRNA_name, current_mRNA_name,
                    filtered_microRNA, filtered_mRNA, lowGroup[current_lncRNA_name].ToString());

                writer.Write("\t" + Rhigh[0] + "\t" + Rhigh[1] + "\t" + Rlow[0] + "\t" + Rlow[1] + "\t");

                //pValue Condition satisfied?
                if (Rhigh[1] > 0.05)
                {
                    writer.Write("p(R_high)");
                    if (Rlow[1] > 0.05)
                        writer.Write(" & p(R_low)");
                    writer.Write(" > 0.05)\t");
                }
                else if (Rlow[1] > 0.05)
                    writer.Write("p(R_low) > 0.05)\t");
                else
                    writer.Write("Satisfied\t");

                //SCC satisfied?
                if (Rlow[0] > -0.4)
                {
                    writer.Write("No (R_low > -0.4)\t");
                }
                else if (Rhigh[0] < -0.4)
                {
                    writer.Write("No (R_high < -0.4)\t");
                }
                else if (Math.Abs(Rlow[0]) - Math.Abs(Rhigh[0]) < 0.4)
                {
                    writer.Write("No (|R_low|-|R_high| < -0.4)\t");
                }
                else
                {
                    //满足条件的triplet
                    //filteredByStep4_triplet.Rows.Add(tripletRow.ItemArray);
                    writer.Write("Satisfied\t");
                }


                //用Fisher transformation变换R为F(R)
                //fhigh=FisherTransformation(Rhigh);
                double Fhigh = FisherTransformation(Rhigh[0]);
                double Flow = FisherTransformation(Rlow[0]);
                writer.Write(Fhigh + "\t" + Flow + "\t");

                //计算rewiring score，调用R，pnorm()函数，该函数给出正态分布随机数小于给定数值的概率。它也被称为“累积分布函数”
                //因为我们的rewiring score是双侧检验，所以得到的p值要*2
                //x <- seq(-10,10,by = .2)
                //y <- pnorm(x, mean = 0, sd = 1)
                int Nhigh = oneForthQuantile, Nlow = oneForthQuantile;
                double miu = Math.Abs(   (Flow-Fhigh)  /  Math.Sqrt(  1.06/(Nhigh-3)+1.06/(Nlow-3)  )   ); //单侧检验，无绝对值，需要low-high>0
                double rewireScore = PValueFromNormalGaussianDistribution(rEngine, miu, true, false);
                writer.Write(rewireScore+"\t");

                //对于当前循环下的triplet，随机生成H-group和L-group，计算rewire score，100次

                //只用找当前lncRNA对应的dataRow
                DataRow dataRow =
                    filtered_lncRNA.Select(filtered_lncRNA.TableName + "='" + current_lncRNA_name + "'")[0];

                //行列转置，但不再需要排序了，直接随机取Nhigh+Nlow个就够了
                DataTable dtSort = new DataTable();
                dtSort.Columns.Add("Sample");
                dtSort.Columns.Add("Expression", typeof(double));
                for (int i = 1;
                    i < lncRNA.Columns.Count;
                    i++) //lncRNA和filtered的两个datatable列结构是一样的，可以直接互用；第0列是lncRNA的Name，不放入表内
                {
                    dtSort.Rows.Add(lncRNA.Columns[i].ColumnName, Convert.ToDouble(dataRow.ItemArray[i]));
                }

                //随机生成H-group和L-group，计算rewire score，100次
                int totalTime = 100;
                int largerThanTrueValueTime = 0;

                for (int randomTime = 0; randomTime < totalTime; randomTime++)
                {
                    //随机一个Nhigh+Nlow的下标列表
                    Random random = new Random();
                    List<int> randomCoordinates = new List<int>();
                    Hashtable duplicatingTest = new Hashtable();
                    while (randomCoordinates.Count < Nhigh + Nlow)
                    {
                        int r = random.Next(0, dtSort.Rows.Count);
                        if (!duplicatingTest.Contains(r))
                        {
                            randomCoordinates.Add(r);
                            duplicatingTest.Add(r, 0);
                        }
                    }

                    //存储random的highGroup和lowGroup的SampleString
                    string randomHighGroupSampleString =
                        dtSort.Rows[randomCoordinates[0]]["Sample"].ToString(); //接sampleString，先放上第0行
                    string randomLowGroupSampleString =
                        dtSort.Rows[randomCoordinates[oneForthQuantile]]["Sample"].ToString();

                    for (int i = 1; i < oneForthQuantile; i++) //接着上面的第0行
                    {
                        randomHighGroupSampleString += "\t" + dtSort.Rows[randomCoordinates[i]]["Sample"];
                        randomLowGroupSampleString += "\t" + dtSort.Rows[randomCoordinates[i + oneForthQuantile]]["Sample"];
                    }
                    
                    double[] randomRhigh = CalculateSCCValueByR(rEngine, current_microRNA_name, current_mRNA_name,
                        filtered_microRNA, filtered_mRNA, randomHighGroupSampleString);

                    double[] randomRlow = CalculateSCCValueByR(rEngine, current_microRNA_name, current_mRNA_name,
                        filtered_microRNA, filtered_mRNA, randomLowGroupSampleString);

                    //计算F
                    double randomFhigh = FisherTransformation(randomRhigh[0]);
                    double randomFlow = FisherTransformation(randomRlow[0]);

                    //计算rewire score
                    double randomMiu = Math.Abs((randomFlow - randomFhigh) / Math.Sqrt(1.06 / (Nhigh - 3) + 1.06 / (Nlow - 3)));
                    double randomRewireScore = PValueFromNormalGaussianDistribution(rEngine, randomMiu, true, false);

                    //判断是否大于real condition的值
                    if (randomRewireScore > rewireScore)
                        largerThanTrueValueTime++;
                }

                double p = 1.0 * largerThanTrueValueTime / totalTime;
                pValues.Add(p);
                writer.WriteLine(p + "\t");

            }
            writer.Close();
            Console.WriteLine("完毕！               ");

            //计算bonferroni-corrected p-adjusted
            Console.Write("计算bonferroni-corrected p值...");
            NumericVector pVector = new NumericVector(rEngine, pValues);
            rEngine.SetSymbol("p", pVector);
            NumericVector p_adjusted = rEngine.Evaluate("p.adjust(p, method = \"bonferroni\", " + pValues.Count + ")")
                .AsNumeric();
            
            //读入并重新写入文件
            List<string> outputDataStrings = new List<string>();

            StreamReader reader = new StreamReader(outputPath4);
            while (!reader.EndOfStream)
                outputDataStrings.Add(reader.ReadLine());
            reader.Close();

            writer = new StreamWriter(outputPath4, false);
            writer.WriteLine(outputDataStrings[0]);
            for (int i = 1; i < outputDataStrings.Count; i++)
            {
                writer.WriteLine(outputDataStrings[i] + p_adjusted[i - 1]);
            }
            writer.Close();
            Console.WriteLine("完毕！");
            //Console.WriteLine("完毕！第4步筛选后的triplet共" + filteredByStep4_triplet.Rows.Count + "行");
            /*
            Console.WriteLine();

            writer = new StreamWriter(@"D:\20190609小肥肥程序数据\step4_filtered_triplet.txt");
            
            foreach (DataRow dataRow in filteredByStep4_triplet.Rows)
            {
                writer.WriteLine(string.Join("\t", dataRow.ItemArray));
            }

            writer.Close();
            */





            Console.WriteLine("按任意键结束程序...");
            Console.ReadLine();
        }
    }
}
