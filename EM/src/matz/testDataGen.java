package matz;

import java.io.*;
import java.util.Date;
import java.util.Random;

public class testDataGen {
	//EMアルゴリズムのテストのために、なんとなく混合正規分布に従っているような感じのテストデータを出力する。
	//Javaの正規乱数(平均0,標準偏差1)を使い、これを各分布ごと適当な倍率scale[k]に拡大、trans[k]平行移動する。
		
	public static int cluster (int clusterMax) {//呼ばれるたびに、設定された数のクラスタのうち、どれかひとつのクラスタを指定して返す
		double rnd = Math.random();
		double minval = Double.MAX_VALUE;
		int index = 0;
		for (int k=0; k<clusterMax; k++) {
			double dist = Math.abs(rnd - (double)k/(double)clusterMax);
			if (dist < minval) {
				minval = dist;
				index = k;
			}
		}
		return index;
	}
	
	public static void main (String args[]) {
		//出力ファイルは改行で区切った1次元ベクトルで、1行に値(時間を想定)のみ書かれる。
		//ファイルは適当な数作れるようにする。
		try {
			int fileNum = 1; //出力ファイルの数
			int clusterMax = 3; //クラスタの数(＝分布の峰の数)
			if (args.length >= 1) { //引数があるなら出力ファイルの数とクラスタの数にそれぞれ対応させる
				fileNum = Integer.parseInt(args[0]);
				if (args.length >= 2) {
					clusterMax = Integer.parseInt(args[1]);
				}
			}
			int size = 20000; 
			
			//実行時にディレクトリがなければ作る
			Date date = new Date();
			String dirNum = Long.toString(date.getTime());
			String path = "data" + dirNum.substring(5);
			File outPath = new File(path);
			if (!outPath.isDirectory()) {
				outPath.mkdir();
			}
			
			int itr = 1;
			int denom = 0;
			while(itr<=clusterMax){ //重み決定用の分母をつくる
				denom += itr;
				itr++;
			}
			
			for (int f = 0; f < fileNum; f++) {
				
				//double[] p = new double[clusterMax]; //重み(所属確率)
				double[] phi = new double[clusterMax];
				double[] trans = new double[clusterMax]; //平行移動幅
				double[] mean = new double[clusterMax];
				double[] scale = new double[clusterMax]; //拡大率
				double[] var = new double[clusterMax]; //分散
				int[] dataNum = new int[clusterMax];
				
				double transMin = 0.0; //平均決定用の下限値
				for (int k=0; k<clusterMax; k++) { //分布のパラメータ初期化
					//p[k] = (clusterMax - k)/denom;
					trans[k] = transMin + 5000*(k+1)/denom + Math.random()*1000 - 500;
					transMin = trans[k];
					scale[k] = Math.random()*2500;
					phi[k] = 0.0;
					var[k] = 0.0;
				}
				
				int[] counter = new int[100];
				for (int c=0; c<100; c++) {
					counter[c] = 0;
				}
				
				//出力
				File outf = new File(path + String.format("/%04d.csv", f));
				OutputStreamWriter osw = new OutputStreamWriter(new FileOutputStream(outf));
				BufferedWriter bw = new BufferedWriter(osw);
				
				double[] data = new double[size];
				int[] cluster = new int[size];
				Random rnd = new Random();
				for (int i=0; i<size; i++) {//データ生成
					int index = cluster(clusterMax);
					cluster[i] = index;
					dataNum[index]++;
					phi[index]++;
					//rnd.setSeed(date.getTime());
					while (true){
						double raw = rnd.nextGaussian();
						data[i] = raw*scale[index] + trans[index];
						if (data[i] >= 0 && data[i] < 10000) {
							break;
						}
					}
					int cIndex = (int)data[i]/100;
					counter[cIndex]++;
					
					mean[index] += data[i];
				}
				
				for (int k=0; k<clusterMax; k++) {
					phi[k] /= size;
					mean[k] /= (double)dataNum[k];
				}
				
				for (int i=0; i<size; i++){
					int index = cluster[i];
					double gap = data[i] - mean[index];
					var[index] += gap*gap;
				}
				
				for (int k=0; k<clusterMax; k++) {
					var[k] /= (double)dataNum[k];
				}
				
				//データ出力
				for (int i=0; i<size; i++) {
					bw.write(String.format("%d,%f", cluster[i], data[i]));
					bw.newLine();
				}
				bw.newLine();
				
				//度数分布表出力
				for (int c=0; c<100; c++) {
					bw.write(String.format("%d,%d", c*100, counter[c]));
					bw.newLine();
				}
				bw.newLine();
				
				//パラメータ出力
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
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}