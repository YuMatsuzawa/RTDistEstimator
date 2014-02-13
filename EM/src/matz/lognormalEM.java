package matz;

/* EMアルゴリズムでモデル推定するクラス。
 * これはLognormalMM（混合対数正規分布）。clusterMaxは本来的にはパラメータなのでここを回して、
 * それも含めて情報量基準か何かで最適モデルを決めないといけないはず。
 * 計算したところ更新式が予定通りGMMと瓜二つだったのでコードも実行時間も瓜二つ。
 * とりあえずRTCountの出力文法に合わせたので、retweetXXXXXXXXXXXXディレクトリ内のXXXXXXXXXXXX.csvファイルをすべて読み込むはず。
 * このLNMMには時間遅れパラメータはまだ導入されていない。
 */

import java.io.*;
import java.util.*;

public class lognormalEM {
	private static int KLEINBERG = 0;
	private static int KMEANS = 1;
	private static int AIC = 0;
	private static int BIC = 1;
	private static int SECONDS_A_WEEK = 60*60*24*7;
	
	private static int clusterMax = 2;
	private static int stateMax = 2;
	private static int maxIter = 10000;
	private static int upperMax = SECONDS_A_WEEK;
	//private static int lowerOffset = 10;
	private static String infroot = "./";
	private static String aiclogFileName = "aic_ln.log";
	private static String biclogFileName = "bic_ln.log";
	private static int clusteringMethod = KLEINBERG; // 0=kleinberg;1=kmeans
	private static int criterion = BIC; //0=AIC;1=BIC
	private static int kleinbergCluster;
	

	
	//private double pseudoZero = 0.001;
	//private double altDelay = 1.0;
	
	public static FilenameFilter regex(String regex){ //正規表現によるファイル名フィルターを返す
		final String regex_ = regex;
		return new FilenameFilter() {
			public boolean accept(File file, String name){
				boolean ret = name.matches(regex_);
				return ret;
			}
		};
	}
	
	//実行時間計測
	private static long start = System.currentTimeMillis();
	private static long split = start;
	private static void stopWatch(String comment) {
		long current = System.currentTimeMillis();
		System.out.println((current - start) + "\t" + (current - split) + "\t" + comment);
		split = current;
	}
	
	public static void main(String args[]) {
		new lognormalEM(args);
		
		try {
			stopWatch("Starting...");
			
			//勝手にデータを探す
			File rootPath = new File(infroot);
			File[] dataDirs = rootPath.listFiles(regex("retweetLog[\\d]+"));
			for (File dir : dataDirs) {
				String[] dataFileNames = dir.list(regex("[\\d]+\\.csv"));
				for (String fileName : dataFileNames) {
					int resFileNum = 1;
					String resFileName = fileName.replaceAll("([\\d]+\\.)(csv)", "$1res\\.ln\\." + resFileNum + "\\.$2");
					File resDir = new File(dir, fileName.replaceAll("([\\d]+\\.)(csv)", "$1res"));
					if(!resDir.isDirectory()) resDir.mkdir();
					File resFile = new File(resDir, resFileName); //結果出力先
					if (clusteringMethod == KMEANS){ //最大クラスタ数を指定して、クラスタ数を増やしながら試行する（K平均法でクラスタ分け）
						for (int cluster=1; cluster <= clusterMax; cluster++){
							while (resFile.exists()){
								//出力ファイルが存在する(既に1回以上処理済み)なら番号をインクリメントする。
								resFileNum++;
								resFileName = fileName.replaceAll("([\\d]+\\.)(csv)", "$1res\\.ln\\." + resFileNum + "\\.$2");
								resFile = new File(resDir, resFileName);
							}
							if (!resFile.exists()) { 
								File dataFile = new File(dir, fileName); //入力元
								new lognormalEM(cluster, dataFile, resFile);
							}
						}
					} else if(clusteringMethod == KLEINBERG) { //バースト検知を使ってクラスタ数を判定する
						while (resFile.exists()){
							//出力ファイルが存在する(既に1回以上処理済み)なら番号をインクリメントする。
							resFileNum++;
							resFileName = fileName.replaceAll("([\\d]+\\.)(csv)", "$1res\\.ln\\." + resFileNum + "\\.$2");
							resFile = new File(resDir, resFileName);
						}
						if (!resFile.exists()){
							File dataFile = new File(dir, fileName); //入力元
							new lognormalEM(clusteringMethod, dataFile, resFile);
						}
					}
					if(criterion == AIC) icLogRename(aiclogFileName, resFile);
					else if(criterion == BIC) icLogRename(biclogFileName, resFile);
				}
			}
			stopWatch("Done.");
		} catch(FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void icLogRename(String iclogFileName, File resFile) {
		File iclogFile = new File(resFile.getParentFile(), iclogFileName);
		int iclogFileNum = 1;
		String icString = iclogFileName.replaceAll("\\.log", "");
		File iclogFileRenameTo = new File(resFile.getParentFile(), resFile.getName().replaceAll("([\\d]+\\.)(res\\.ln\\.[\\d]+\\.)(csv)", "$1" + icString +"\\." + iclogFileNum + "\\.$3"));
		while (iclogFileRenameTo.exists()){
			//出力ファイルが存在する(既に1回以上処理済み)なら番号をインクリメントする。
			iclogFileNum++;
			iclogFileRenameTo = new File(resFile.getParentFile(), resFile.getName().replaceAll("([\\d]+\\.)(res\\.ln\\.[\\d]+\\.)(csv)", "$1" + icString +"\\." + iclogFileNum + "\\.$3"));
		}
		if(iclogFile.length() > 0) {
			iclogFile.renameTo(iclogFileRenameTo);
			stopWatch(iclogFileRenameTo.toString());
		}
		
	}
	
	public lognormalEM(String[] args) { //初期化コンストラクタ
		try {
			/* 引数の順序は固定とする．
			 * デフォルトから変更しない後続の引数は省略してもいいが，
			 * 前方の引数を省略して後方の引数だけを入力することはできない．
			 * 後方の引数だけを変更したい場合は前方の引数のデフォルト値を表記する．
			 * デフォルトは"./ 2 10000 0 1"
			 */
			if (args.length >= 1) {
				infroot = args[0];
				if (args.length >= 2) {
					clusterMax = Integer.parseInt(args[1]);
					stateMax = Integer.parseInt(args[1]);
					if(stateMax < 2) stateMax = 2;
					if (args.length >= 3) {
						maxIter = Integer.parseInt(args[2]);
						if (args.length >= 4) {
							upperMax = Integer.parseInt(args[3]);
							if (args.length >= 5) {
								clusteringMethod = Integer.parseInt(args[4]);
								if (args.length >= 6) {
									criterion = Integer.parseInt(args[5]);
								}
							}
						}
					}
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
	
	public lognormalEM(int cluster, File dataFile, File resFile) throws FileNotFoundException, IOException { //計算コンストラクタ。見ての通り入力データファイルと出力先ファイルを指定する。
		List<Double> dataList = new ArrayList<Double>();
		
		stopWatch("Reading: " + dataFile.toString());
		/* RTログの方は1行目にオリジナルTweetのURLが書かれているので、データは2行目から入る。
		 * よってそのオフセットを指定しておくこと。
		 */

		InputStreamReader isr = new InputStreamReader(new FileInputStream(dataFile));
		BufferedReader br = new BufferedReader(isr);
		String head = br.readLine(); //先頭行
		
		String line = new String();
		while((line = br.readLine()) != null) {//データをリストに読み込み。null（ファイル末尾）でブレイクする。
			try {
				StringTokenizer st = new StringTokenizer(line, ",");
				dataList.add(Double.parseDouble(st.nextToken()));
			} catch(NumberFormatException e) { //Doubleにパースできない変な値が見つかったらリターンしておく。
				e.printStackTrace();
				br.close();
				isr.close();
				stopWatch("Unexpected Data Type. Aborting.");
				return;
			}
		}
		br.close();
		isr.close();
		
		int size = 0;
		for(int i=0; i<dataList.size(); i++) if(dataList.get(i) <= upperMax) size++;
		if (size < 100) { //1週間以内のRT数が100未満なら抜ける。若干特殊すぎる
			stopWatch("Disqualified. Aborting.");
			return;
		}
		//double[] data = new double[dataList.size()];
		double[] data = new double[size];
		int dataMin = Integer.MAX_VALUE;
		int dataMax = 0;
		int dataIndex = 0;
		for(int i=0; i<dataList.size(); i++){ //ここでリストを配列に移し替え、以降はこの配列からデータを参照する。doubleがプリミティブ型なのでちょうどいいメソッドがないらしい
			if (dataList.get(i) <= upperMax) {//使用するデータを極限で絞る
				data[dataIndex] = dataList.get(i);
				//data[dataIndex] = dataList.get(i) + lowerOffset;
				if (data[dataIndex] < dataMin) dataMin = (int)data[dataIndex];
				if (data[dataIndex] > dataMax) dataMax = (int)data[dataIndex];
				dataIndex++;
			}
		}
		dataList.clear();
		
		//初期所属クラスタを決定
		int[] initCluster = new int[size];
		if(cluster > 0) initCluster = kmeans(cluster, data, dataMin, dataMax, maxIter);
		else if(cluster == 0) {
			initCluster = kleinberg(cluster, data, dataMin, dataMax, maxIter);
			cluster = kleinbergCluster;
		}
		stopWatch("Processing EM for " + cluster + " cluster..");
		
		//パラメータの初期化
		double[] weight = new double[cluster];
		double[] mean = new double[cluster];
		double[] var = new double[cluster];
		
		for (int k=0; k<cluster; k++) {
			weight[k] = 0.0;
			mean[k] = 0.0;
			var[k] = 0.0;
			
		}
		int[] dataNum = new int[cluster];
		for(int i=0; i<size; i++){
			int index = initCluster[i];
			dataNum[index]++;
			weight[index]++;
			/*if(data[i] >= 1) mean[index] += Math.log(data[i]);   //!!!
			else mean[index] += Math.log(altDelay);*/
			mean[index] += Math.log(data[i]);
		}
		for(int k=0; k<cluster; k++) {
			if(dataNum[k]==0) continue;
			weight[k] /= (double)size;
			mean[k] /= (double)dataNum[k];
		}
		for(int i=0; i<size; i++) {
			int index = initCluster[i];
			double gap = Math.log(data[i]) - mean[index];
			/*if(data[i] >= 1) gap = Math.log(data[i]) - mean[index];
			else gap = Math.log(altDelay) - mean[index];*/
			var[index] += gap*gap;
		}
		for(int k=0; k<cluster; k++) {
			var[k] /= (double)dataNum[k];
		}
		
		//イテレータ
		double[][] eAff = new double[size][cluster];
		final double cnst = Math.sqrt(2.0*Math.PI);
		ArrayList<Double> likelihoodSeries = new ArrayList<Double>();
		for(int itr=0; itr<maxIter; itr++){
			
			//Eステップ
			double[] denom = new double[size];
			for(int i=0; i<size; i++){
				denom[i] = 0.0;
				for(int k=0; k<cluster; k++){
					try {
						if(var[k]==0) throw new Exception();
						/*else if(data[i] == 0){
							double gap = Math.log(altDelay) - mean[k];
							eAff[i][k] = weight[k] * Math.exp(-0.5*gap*gap/var[k])/(cnst*Math.sqrt(var[k]))/altDelay;
						}*/else {
							double gap = Math.log(data[i]) - mean[k];
							eAff[i][k] = weight[k] * Math.exp(-0.5*gap*gap/var[k])/(cnst*Math.sqrt(var[k]))/data[i];
						}
					} catch (Exception e) { //不適切な分散を持つクラスタをcatch
						e.printStackTrace();
						stopWatch("Some cluster doesn't have sufficient variance.");
						return;
					}
					denom[i] += eAff[i][k];
				}
				for(int k=0; k<cluster; k++){
					eAff[i][k] /= denom[i];
				}
			}
			
			//Mステップ
			double[] sumAff = new double[cluster];
			for(int k=0; k<cluster; k++){
				weight[k] = 0.0;
				sumAff[k] = 0.0;
				mean[k] = 0.0;
				var[k] = 0.0;
				
				for(int i=0; i<size; i++){
					sumAff[k] += eAff[i][k];
					/*if(data[i]!=0) mean[k] += Math.log(data[i])*eAff[i][k];
					else mean[k] += Math.log(altDelay)*eAff[i][k];*/
					mean[k] += Math.log(data[i])*eAff[i][k];
				}
				weight[k] = sumAff[k]/(double)size;
				mean[k] /= sumAff[k];
				
				for(int i=0; i<size; i++){
					double gap = Math.log(data[i]) - mean[k];
					var[k] += gap*gap*eAff[i][k];
					/*if(data[i]!=0) {
						double gap = Math.log(data[i]) - mean[k];
						var[k] += gap*gap*eAff[i][k];
					}else{
						double gap = Math.log(altDelay) - mean[k];
						var[k] += gap*gap*eAff[i][k];
					}*/
				}
				
				var[k] /= sumAff[k];
				
				try{
					if(Double.isNaN(weight[k]) || Double.isNaN(var[k])){
						System.out.println(k);
						throw new Exception();
					}
				}catch(Exception e){
					e.printStackTrace();
					return;
				}
			}
			
			//尤度
			double maxLikelihood = 0.0;
			for(int i=0; i<size; i++){
				double dataLikelihood = 0.0;
				for (int k=0; k<cluster; k++){
					/*double delay = data[i];
					if (delay < 1) delay = altDelay;
					double gap = Math.log(delay) - mean[k];
					double denomK = Math.sqrt(2 * Math.PI * var[k]) * delay;*/
					double gap = Math.log(data[i]) - mean[k];
					double denomK = Math.sqrt(2 * Math.PI * var[k]) * data[i];
					dataLikelihood += weight[k] * Math.exp(-0.5 * gap * gap / var[k]) / denomK;
				}
				maxLikelihood += Math.log(dataLikelihood);
			}
			if(itr % 10 == 0){
				likelihoodSeries.add(maxLikelihood);
				//stopWatch("\tLikelihood:\t" + maxLikelihood);
			}
		}
		
		
		stopWatch("Writing: " + resFile.toString());
		//結果出力
		OutputStreamWriter osw = new OutputStreamWriter(new FileOutputStream(resFile));
		BufferedWriter bw = new BufferedWriter(osw);
		bw.write(head);
		bw.newLine();
		
		//近似曲線プロット生成・出力
		int range = dataMax - dataMin;
		int pitch = 10*60;
		int resol = range/pitch;
		for (int r=0; r<resol; r++) {
			double plot = 0.0;
			double subPlot = 0.0;
			int time = dataMin + pitch*r;
			bw.write(String.format("%d", time));
			for (int k=0; k<cluster; k++) {
				subPlot = weight[k] * lognormal(time, mean[k], var[k]);
				bw.write(String.format(",%f", subPlot));
				plot += subPlot;
			}
			bw.write(String.format(",%f", plot));
			bw.newLine();
		}
		bw.newLine();
		
		//クラスタ成分出力
		for(int k=0; k<cluster; k++) {
			bw.write(String.format("#Cluster No.%d", k));
			bw.newLine();
			bw.write(String.format("weight: %f", weight[k]));
			bw.newLine();
			bw.write(String.format("mean: %f", mean[k]));
			bw.newLine();
			bw.write(String.format("variance: %f", var[k]));
			bw.newLine();
			bw.newLine();
		}
		
		//尤度推移出力
		for(int i=0; i<likelihoodSeries.size(); i++){
			bw.write(likelihoodSeries.get(i).toString());
			bw.newLine();
		}
		bw.close();
		osw.close();
		
		//criterion出力
		if (criterion == AIC) AIC(cluster, size, data, mean, var, weight, resFile);
		else if (criterion == BIC) BIC(cluster, size, data, mean, var, weight, resFile);
	}
	
	private void BIC(int cluster, int size, double[] data, double[] mean,
			double[] var, double[] weight, File resFile) {//BIC出力
		try {
			File biclogFile = new File(resFile.getParentFile(), biclogFileName);
			OutputStreamWriter biclogOsw = new OutputStreamWriter(new FileOutputStream(biclogFile, true));
			BufferedWriter biclogBw = new BufferedWriter(biclogOsw);
			double maxLikelihood = 0.0;
			for(int i=0; i<size; i++){
				double dataLikelihood = 0.0;
				for (int k=0; k<cluster; k++){
					/*double delay = data[i];
					if (delay < 1) delay = altDelay;
					double gap = Math.log(delay) - mean[k];
					double denom = Math.sqrt(2 * Math.PI * var[k]) * delay;*/
					double gap = Math.log(data[i]) - mean[k];
					double denom = Math.sqrt(2 * Math.PI * var[k]) * data[i];
					dataLikelihood += weight[k] * Math.exp(-0.5 * gap * gap / var[k]) / denom;
				}
				maxLikelihood += Math.log(dataLikelihood);
			}
			int numOfParam = 3; //クラスタあたりのパラメータ数
			double BIC = -2 * maxLikelihood + Math.log(size) * (numOfParam * cluster + cluster - 1);
			
			if (biclogFile.length() == 0){
				stopWatch("Using BIC...");
				biclogBw.write("maxIter = ," + maxIter + ", upperMax = ," + upperMax);
				biclogBw.newLine();
			}
			biclogBw.write(cluster + ", " + BIC + ", " + resFile.toString());
			biclogBw.newLine();
			
			biclogBw.close();
			biclogOsw.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return;
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}

	private void AIC(int cluster, int size, double[] data, double[] mean,
			double[] var, double[] weight, File resFile) { //AIC出力
		try {
			File aiclogFile = new File(resFile.getParentFile(), aiclogFileName);
			OutputStreamWriter aiclogOsw = new OutputStreamWriter(new FileOutputStream(aiclogFile, true));
			BufferedWriter aiclogBw = new BufferedWriter(aiclogOsw);
			double maxLikelihood = 0.0;
			for(int i=0; i<size; i++){
				double dataLikelihood = 0.0;
				for (int k=0; k<cluster; k++){
					/*double delay = data[i];
					if (delay < 1) delay = altDelay;
					double gap = Math.log(delay) - mean[k];
					double denom = Math.sqrt(2 * Math.PI * var[k]) * delay;*/
					double gap = Math.log(data[i]) - mean[k];
					double denom = Math.sqrt(2 * Math.PI * var[k]) * data[i];
					dataLikelihood += weight[k] * Math.exp(-0.5 * gap * gap / var[k]) / denom;
				}
				maxLikelihood += Math.log(dataLikelihood);
			}
			int numOfParam = 3; //クラスタあたりのパラメータ数
			double AIC = -2 * maxLikelihood + 2 * (numOfParam * cluster + cluster - 1);
			
			if (aiclogFile.length() == 0){
				stopWatch("Using AIC...");
				aiclogBw.write("maxIter = ," + maxIter + ", upperMax = ," + upperMax);
				aiclogBw.newLine();
			}
			aiclogBw.write(cluster + ", " + AIC + ", " + resFile.toString());
			aiclogBw.newLine();
			
			aiclogBw.close();
			aiclogOsw.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return;
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
	
/*	private double delayedLognormal(int time, double tau, double mean, double var) {
		double delay = time - tau;
		if(delay <= 0) return 0.0;
		else if(delay < 1) delay = 1;
		double denom = Math.sqrt(2*Math.PI*var) * delay;
		double gap = Math.log(delay) - mean;
		return Math.exp(-0.5*gap*gap/var)/denom;
	}*/
	
	private double lognormal(int time, double mean, double var) {
		double denom = Math.sqrt(2*Math.PI*var) * time;
		double gap = Math.log(time) - mean;
		return Math.exp(-0.5*gap*gap/var)/denom;
	}

/*	public double gaussian(int time, double mean, double var) {
		double denom = Math.sqrt(2*Math.PI*var);
		double gap = time - mean;
		return Math.exp(-0.5*gap*gap/var)/denom;
	}*/

	public int[] kmeans(int cluster, double[] data, double dataMin, double dataMax, int max_itr) {//k平均法。各データ点が最初に属するクラスタのインデックスを値として持つ配列を返す。
		stopWatch("Clustering by K-means method..");
		
		int size = data.length;
		int[] initCluster = new int[size];
		
		//クラスタ中心
		double[] centers = new double[cluster];
		for (int k=0; k<cluster; k++) {
			centers[k] = dataMin + k*(dataMax - dataMin)/cluster; //レンジ内に等間隔に配置
			//centers[k] = dataMin + (dataMax - dataMin)*Math.random(); //レンジ内でランダムに配置
		}
		
		//イテレータ
		for (int itr=0; itr<max_itr; itr++) {
			int cutoff = size/100;
			
			int[] dataNum2 = new int[cluster]; 
			/* 最初にランダム生成された際、位置が悪くて（集団から遠くて）メンバが一人以下しか割り振られなかったクラスタがそのままになっていると、後にバグを引き起こす（そのクラスタだけ分散が0になってしまうので）。
			 * これを避けるためにもしメンバが一人以下のクラスタが見つかったらそのクラスタだけ中心を恣意的に調整する。
			 */
			//最も中心の近いクラスタにデータを配属
			for(int i=0; i<size; i++) {
				double minval = Double.MAX_VALUE;
				int index = 0;
				for(int k=0; k<cluster; k++){
					double dist = Math.abs(data[i] - centers[k]);
					if(dist < minval) {
						minval = dist;
						index = k;
					}
					initCluster[i] = index;
				}
				dataNum2[index]++;
			}
			
			//中心再計算
			int[] dataNum = new int[cluster];
			for(int k=0; k<cluster; k++){
				centers[k] = 0.0;
			}
			for(int i=0; i<size; i++) {
				int index = initCluster[i];
				centers[index] += data[i];
				dataNum[index]++;
			}
			for(int k=0; k<cluster; k++) {
				if(dataNum[k]<=cutoff) continue; //落ちこぼれは後で救済する。落ちこぼれ閾値cutoffを導入
				centers[k] /= dataNum[k];
			}
			for(int k=0; k<cluster; k++) {//落ちこぼれ救済
				double delta = 0.001 * (dataMax - dataMin);
				if(dataNum[k]<=cutoff) {
					int l = k;
					while (l==k || dataNum[l] <= cutoff) l = (int)(Math.random()*cluster);
					centers[k] = centers[l] + 0.5*delta;
					centers[l] = centers[l] - 0.5*delta; //既にメンバを持っているクラスタをランダムに取ってきて、そのクラスタの中心を中心に幅deltaだけ離れた位置に中心を取らせてもらうことで救済する
				}
			}
		}
		return initCluster;
	}
	
	public int[] kleinberg(int cluster, double[] data, double dataMin, double dataMax, int max_itr){ //Kleinbergのバースト検知アルゴリズムによるクラスタリング
		stopWatch("Clustering by Kleinberg Algorithm..");
		
		int size = data.length;
		
		//frequency map
		double range = dataMax - dataMin;
		//int pitch = 10; //頻度分布の区間（分）
		int pitch = 10*60; //頻度分布の区間（秒）
		int arraySize = 1 + (int)range / pitch;
		int[] dataFreq = new int[arraySize];
		int[] state = new int[arraySize];
		for(int a=0; a<arraySize; a++) {
			dataFreq[a] = 0;
			state[a] = 0;
		}
		for(int i=0; i<size; i++){
			int index = (int)(((int)data[i] - (int)dataMin) / pitch);
			dataFreq[index]++;
		}
		int maxFreq = 0;
		int minFreq = Integer.MAX_VALUE;
		double avgFreq = 0.0;
		for(int a=0; a<arraySize; a++){
			avgFreq += dataFreq[a];
			if(dataFreq[a] > maxFreq) maxFreq = dataFreq[a];
			if(dataFreq[a] < minFreq) minFreq = dataFreq[a];
		}
		avgFreq /= arraySize;
		
		//kleinberg
		int[] initCluster = new int[size];

		//double level = 1.0;
		double gamma = 1.0;
		double numer = avgFreq - minFreq;
		double denom = maxFreq - minFreq;
		int denomInt = 10;
		double level = (denom / numer - 1) / stateMax;
		//double level = (denom / numer) / stateMax;
		//if(level > denom/numer) level = denom / numer;
		double normal = numer / denom;
		double[] stateCost = new double[stateMax];
		for (int s=0; s<stateMax; s++) {
			stateCost[s] = normal * (1 + level * s);
			//stateCost[s] = normal * (level * s);
		}
		
		double sequentialCost = 0.0;		
		for(int a=0; a<arraySize; a++){
			int freqInt = (int) (denomInt * (dataFreq[a] - minFreq) / denom);
			//状態コスト
			double[] candidateCost = new double[stateMax];
			for(int s=0; s<stateMax; s++){
				candidateCost[s] = - Math.log(Math.pow(stateCost[s], freqInt) * Math.pow(1 - stateCost[s], denomInt - freqInt) * combination(denomInt, freqInt))
						+ sequentialCost;
				if(a>1 && state[a-1] < s) candidateCost[s] += (s - state[a-1]) * gamma; //状態遷移コスト
			}

			double lowestCost = Double.MAX_VALUE;
			int lowestCostState = 0;
			for(int s=0; s<stateMax; s++){
				if(candidateCost[s] < lowestCost) {
					lowestCost = candidateCost[s];
					lowestCostState = s;
				}
			}
			sequentialCost = lowestCost;
			state[a] = lowestCostState;
		}
		
		ArrayList<Integer> burstPoints = new ArrayList<Integer>();
		burstPoints.add(0);
		boolean bursting = false;
		if (state[0] > 0) bursting = true;
		for(int a=1; a<arraySize; a++){
			if(state[a] > state[a-1] && !bursting) {
				if(!burstPoints.contains(a-1)) burstPoints.add(a-1);
				bursting = true;
			}else if(state[a] < state[a-1] && bursting){
				bursting = false;
			}
		}
		
		for(int k=0; k<burstPoints.size(); k++){
			for(int i=0; i<size; i++){
				if(data[i] - dataMin >= burstPoints.get(k) * pitch) initCluster[i] = k;
			}
		}
		
		kleinbergCluster = burstPoints.size();
		
		return initCluster;
	}
	
	public int combination(int n, int k) {
		double numer = 1.0;
		double denom = 1.0;
		int start;
		int end;
		if(k > n/2){
			start = k+1;
			end = n-k;
		}else{
			start = n-k+1;
			end = k;
		}
		for (int i=start; i<=n; i++) numer *= i;
		for (int i=1; i<=end; i++) denom *= i;
		return (int) (numer/denom);
		
	}
}



