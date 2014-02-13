package matz;

/* EMアルゴリズムでモデル推定するクラス。
 * これはGMM（混合正規分布）。clusterMaxは本来的にはパラメータなのでここを回して、
 * それも含めて情報量基準か何かで最適モデルを決めないといけないはず。
 * このクラスは実行時間はそれほど大きくない。ただし、clusterMaxを大きくすると鬼のような実行時間になると思われる
 * とりあえずRTCountの出力文法に合わせたので、retweetXXXXXXXXXXXXディレクトリ内のXXXXXXXXXXXX.csvファイルをすべて読み込むはず。
 */

import java.io.*;
import java.util.*;

public class EM {
	private static int clusterMax = 3;
	private static int maxIter = 30;
	private static String infroot = "./";
	
	public static FilenameFilter regex(String regex){ //正規表現によるファイル名フィルターを返す
		final String regex_ = regex;
		return new FilenameFilter() {
			public boolean accept(File file, String name){
				boolean ret = name.matches(regex_);
				return ret;
			}
		};
	}
	
	
	public static void main(String args[]) {
		
		try { //コンストラクタ的なことは先にここで行う。
			if (args.length >= 1) {
				infroot = args[0];
				if (args.length >= 2) {
					clusterMax = Integer.parseInt(args[1]);
					if (args.length >= 3) {
						maxIter = Integer.parseInt(args[2]);
					}
				}
			}
		} catch (Exception e) {
			System.out.println("Usage: matz.EM[ <ROOT_DIR>[ <#_OF_CLUSTER>[ <#_OF_ITERATION>]]]");
			System.out.println("Note: Arguments must follow the order above.");
			e.printStackTrace();
			return;
		}
		
		try {
			//if (args.length >= 1) clusterMax = Integer.parseInt(args[0]); //この初期設定は上のtry文節に移動した
			//if (args.length >= 2) maxIter = Integer.parseInt(args[1]);
			
			//勝手にデータを探す
			File rootPath = new File(infroot);
			File[] dataDirs = rootPath.listFiles(regex("retweetLog[\\d]+"));
			for (File dir : dataDirs) {
				String[] dataFileNames = dir.list(regex("[\\d]+\\.csv"));
				for (String fileName : dataFileNames) {
					File resFile = new File(dir, fileName.replaceAll("([\\d]+\\.)(csv)", "$1res.$2")); //結果出力先
					if (!resFile.exists()) { 
						//出力ファイルが存在する(既に処理済みの)データは無視する
						File dataFile = new File(dir, fileName); //入力元
						new EM(dataFile, resFile);
					}
				}
			}
		} catch(FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
		
	public EM(File dataFile, File resFile) throws FileNotFoundException, IOException { //計算コンストラクタ。見ての通り入力データファイルと出力先ファイルを指定する。
		List<Double> dataList = new ArrayList<Double>();
		
		
		/* RTログの方は1行目にオリジナルTweetのURLが書かれているので、データは2行目から入る。
		 * よってそのオフセットを指定しておくこと。
		 */

		InputStreamReader isr = new InputStreamReader(new FileInputStream(dataFile));
		BufferedReader br = new BufferedReader(isr);
		String head = br.readLine(); //先頭行
		
		String line = new String();
		while((line = br.readLine()) != null) {//データをリストに読み込み。null（ファイル末尾）でブレイクする。
			//String[] splitedLine = line.split(","); //スプリッターはテストデータでは使っていたがRTログでは使わない
			//dataList.add(Double.parseDouble(splitedLine[1]));
			try {
				dataList.add(Double.parseDouble(line));
			} catch(NumberFormatException e) { //Doubleにパースできない変な値が見つかったらリターンしておく。
				e.printStackTrace();
				br.close();
				isr.close();
				return;
			}
		}
		br.close();
		isr.close();
			
		double[] data = new double[dataList.size()];
		int dataMin = (int)Double.MAX_VALUE; //あとで近似曲線を表示させる時のために最大値と最小値を取っておく。小数点以下はいらない。
		int dataMax = 0;
		for(int i=0; i<dataList.size(); i++){ //ここでリストを配列に移し替え、以降はこの配列からデータを参照する。doubleがプリミティブ型なのでちょうどいいメソッドがないらしい
			data[i] = dataList.get(i);
			if (data[i] < dataMin) dataMin = (int)data[i];
			if (data[i] > dataMax) dataMax = (int)data[i];
		}
		int size = data.length;
		
		//k平均法でデータの初期所属クラスタを決定
		int[] initCluster = kmeans(data, dataMin, dataMax, maxIter);
		
		//Explicitパラメータの初期化
		double[] phi = new double[clusterMax];
		double[] mean = new double[clusterMax];
		double[] var = new double[clusterMax];
		for (int k=0; k<clusterMax; k++) {
			phi[k] = 0.0;
			mean[k] = 0.0;
			var[k] = 0.0;
		}
		int[] dataNum = new int[clusterMax];
		for(int i=0; i<size; i++){
			int index = initCluster[i];
			dataNum[index]++;
			phi[index]++;
			mean[index] += data[i];
		}
		for(int k=0; k<clusterMax; k++) {
			if(dataNum[k]==0) continue;
			phi[k] /= (double)size;
			mean[k] /= (double)dataNum[k];
		}
		for(int i=0; i<size; i++) {
			int index = initCluster[i];
			double gap = data[i] - mean[index];
			var[index] += gap*gap;
		}
		for(int k=0; k<clusterMax; k++) {
			var[k] /= (double)dataNum[k];
		}
		
		//イテレータ
		double[][] eAff =new double[size][clusterMax];
		double[] sqrtVar = new double[clusterMax];
		final double cnst = Math.sqrt(2.0*Math.PI);
		for(int itr=0; itr<maxIter; itr++){
			
			for(int k=0; k<clusterMax; k++){
				sqrtVar[k] = Math.sqrt(var[k]);
			}
			
			//Eステップ
			double[] denom = new double[size];
			for(int i=0; i<size; i++){
				denom[i] = 0.0;
				for(int k=0; k<clusterMax; k++){
					double gap = data[i] - mean[k];
					eAff[i][k] = phi[k] * Math.exp(-0.5*gap*gap/var[k])/(cnst*sqrtVar[k]);
					denom[i] += eAff[i][k];
				}
				for(int k=0; k<clusterMax; k++){
					eAff[i][k] /= denom[i];
				}
			}
			
			//Mステップ
			double[] sumAff = new double[clusterMax];
			for(int k=0; k<clusterMax; k++){
				phi[k] = 0.0;
				sumAff[k] = 0.0;
				mean[k] = 0.0;
				var[k] = 0.0;
				
				for(int i=0; i<size; i++){
					sumAff[k] += eAff[i][k];
					mean[k] += data[i]*eAff[i][k];
				}
				phi[k] = sumAff[k]/(double)size;
				mean[k] /= sumAff[k];
				
				for(int i=0; i<size; i++){
					double gap = data[i] - mean[k];
					var[k] += gap*gap*eAff[i][k];
				}
				
				var[k] /= sumAff[k];
			}
		}
		
		
		//結果出力
		OutputStreamWriter osw = new OutputStreamWriter(new FileOutputStream(resFile));
		BufferedWriter bw = new BufferedWriter(osw);
		bw.write(head);
		bw.newLine();
		
		//近似曲線プロット生成・出力
		int range = dataMax - dataMin;
		int pitch = 10;
		int resol = range/pitch;
		for (int r=0; r<resol; r++) {
			double plot = 0.0;
			double subPlot = 0.0;
			int time = dataMin + pitch*r;
			bw.write(String.format("%d", time));
			for (int k=0; k<clusterMax; k++) {
				subPlot = phi[k]*gaussian(time,mean[k],var[k]);
				bw.write(String.format(",%f", subPlot));
				plot += subPlot;
			}
			bw.write(String.format(",%f", plot));
			bw.newLine();
		}
		bw.newLine();
		
		//クラスタ成分出力
		for(int k=0; k<clusterMax; k++) {
			bw.write(String.format("#Cluster No.%d", k));
			bw.newLine();
			bw.write(String.format("weight: %f", phi[k]));
			bw.newLine();
			bw.write(String.format("mean: %f", mean[k]));
			bw.newLine();
			bw.write(String.format("variance: %f", var[k]));
			bw.newLine();
			bw.newLine();
		}
		bw.close();
		osw.close();
		
		
	}
	
	public double gaussian(int time, double mean, double var) {
		double denom = Math.sqrt(2*Math.PI*var);
		double gap = time - mean;
		return Math.exp(-0.5*gap*gap/var)/denom;
	}

	public int[] kmeans(double[] data, double dataMin, double dataMax, int max_itr) {//k平均法。各データ点が最初に属するクラスタのインデックスを値として持つ配列を返す。
		int size = data.length;
		int[] cluster = new int[size];
		
		//クラスタ中心
		double[] centers = new double[clusterMax];
		for (int k=0; k<clusterMax; k++) {
			centers[k] = dataMin + (dataMax - dataMin)*Math.random();
		}
		
		//イテレータ
		for (int itr=0; itr<max_itr; itr++) {
			//最も中心の近いクラスタにデータを配属
			for(int i=0; i<size; i++) {
				double minval = Double.MAX_VALUE;
				int index = 0;
				for(int k=0; k<clusterMax; k++){
					double dist = Math.abs(data[i] - centers[k]);
					if(dist < minval) {
						minval = dist;
						index = k;
					}
					cluster[i] = index;
				}
			}
			
			//中心再計算
			int[] dataNum = new int[clusterMax];
			for(int k=0; k<clusterMax; k++){
				centers[k] = 0.0;
			}
			for(int i=0; i<size; i++) {
				int index = cluster[i];
				centers[index] += data[i];
				dataNum[index]++;
			}
			for(int k=0; k<clusterMax; k++) {
				if(dataNum[k]==0) continue;
				centers[k] /= dataNum[k];
			}
		}
		return cluster;
	}
}