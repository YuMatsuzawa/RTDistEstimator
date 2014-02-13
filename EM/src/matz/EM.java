package matz;

/* EM�A���S���Y���Ń��f�����肷��N���X�B
 * �����GMM�i�������K���z�j�BclusterMax�͖{���I�ɂ̓p�����[�^�Ȃ̂ł������񂵂āA
 * ������܂߂ď��ʊ�������ōœK���f�������߂Ȃ��Ƃ����Ȃ��͂��B
 * ���̃N���X�͎��s���Ԃ͂���قǑ傫���Ȃ��B�������AclusterMax��傫������ƋS�̂悤�Ȏ��s���ԂɂȂ�Ǝv����
 * �Ƃ肠����RTCount�̏o�͕��@�ɍ��킹���̂ŁAretweetXXXXXXXXXXXX�f�B���N�g������XXXXXXXXXXXX.csv�t�@�C�������ׂēǂݍ��ނ͂��B
 */

import java.io.*;
import java.util.*;

public class EM {
	private static int clusterMax = 3;
	private static int maxIter = 30;
	private static String infroot = "./";
	
	public static FilenameFilter regex(String regex){ //���K�\���ɂ��t�@�C�����t�B���^�[��Ԃ�
		final String regex_ = regex;
		return new FilenameFilter() {
			public boolean accept(File file, String name){
				boolean ret = name.matches(regex_);
				return ret;
			}
		};
	}
	
	
	public static void main(String args[]) {
		
		try { //�R���X�g���N�^�I�Ȃ��Ƃ͐�ɂ����ōs���B
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
			//if (args.length >= 1) clusterMax = Integer.parseInt(args[0]); //���̏����ݒ�͏��try���߂Ɉړ�����
			//if (args.length >= 2) maxIter = Integer.parseInt(args[1]);
			
			//����Ƀf�[�^��T��
			File rootPath = new File(infroot);
			File[] dataDirs = rootPath.listFiles(regex("retweetLog[\\d]+"));
			for (File dir : dataDirs) {
				String[] dataFileNames = dir.list(regex("[\\d]+\\.csv"));
				for (String fileName : dataFileNames) {
					File resFile = new File(dir, fileName.replaceAll("([\\d]+\\.)(csv)", "$1res.$2")); //���ʏo�͐�
					if (!resFile.exists()) { 
						//�o�̓t�@�C�������݂���(���ɏ����ς݂�)�f�[�^�͖�������
						File dataFile = new File(dir, fileName); //���͌�
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
		
	public EM(File dataFile, File resFile) throws FileNotFoundException, IOException { //�v�Z�R���X�g���N�^�B���Ă̒ʂ���̓f�[�^�t�@�C���Əo�͐�t�@�C�����w�肷��B
		List<Double> dataList = new ArrayList<Double>();
		
		
		/* RT���O�̕���1�s�ڂɃI���W�i��Tweet��URL��������Ă���̂ŁA�f�[�^��2�s�ڂ������B
		 * ����Ă��̃I�t�Z�b�g���w�肵�Ă������ƁB
		 */

		InputStreamReader isr = new InputStreamReader(new FileInputStream(dataFile));
		BufferedReader br = new BufferedReader(isr);
		String head = br.readLine(); //�擪�s
		
		String line = new String();
		while((line = br.readLine()) != null) {//�f�[�^�����X�g�ɓǂݍ��݁Bnull�i�t�@�C�������j�Ńu���C�N����B
			//String[] splitedLine = line.split(","); //�X�v���b�^�[�̓e�X�g�f�[�^�ł͎g���Ă�����RT���O�ł͎g��Ȃ�
			//dataList.add(Double.parseDouble(splitedLine[1]));
			try {
				dataList.add(Double.parseDouble(line));
			} catch(NumberFormatException e) { //Double�Ƀp�[�X�ł��Ȃ��ςȒl�����������烊�^�[�����Ă����B
				e.printStackTrace();
				br.close();
				isr.close();
				return;
			}
		}
		br.close();
		isr.close();
			
		double[] data = new double[dataList.size()];
		int dataMin = (int)Double.MAX_VALUE; //���Ƃŋߎ��Ȑ���\�������鎞�̂��߂ɍő�l�ƍŏ��l������Ă����B�����_�ȉ��͂���Ȃ��B
		int dataMax = 0;
		for(int i=0; i<dataList.size(); i++){ //�����Ń��X�g��z��Ɉڂ��ւ��A�ȍ~�͂��̔z�񂩂�f�[�^���Q�Ƃ���Bdouble���v���~�e�B�u�^�Ȃ̂ł��傤�ǂ������\�b�h���Ȃ��炵��
			data[i] = dataList.get(i);
			if (data[i] < dataMin) dataMin = (int)data[i];
			if (data[i] > dataMax) dataMax = (int)data[i];
		}
		int size = data.length;
		
		//k���ϖ@�Ńf�[�^�̏��������N���X�^������
		int[] initCluster = kmeans(data, dataMin, dataMax, maxIter);
		
		//Explicit�p�����[�^�̏�����
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
		
		//�C�e���[�^
		double[][] eAff =new double[size][clusterMax];
		double[] sqrtVar = new double[clusterMax];
		final double cnst = Math.sqrt(2.0*Math.PI);
		for(int itr=0; itr<maxIter; itr++){
			
			for(int k=0; k<clusterMax; k++){
				sqrtVar[k] = Math.sqrt(var[k]);
			}
			
			//E�X�e�b�v
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
			
			//M�X�e�b�v
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
		
		
		//���ʏo��
		OutputStreamWriter osw = new OutputStreamWriter(new FileOutputStream(resFile));
		BufferedWriter bw = new BufferedWriter(osw);
		bw.write(head);
		bw.newLine();
		
		//�ߎ��Ȑ��v���b�g�����E�o��
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
		
		//�N���X�^�����o��
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

	public int[] kmeans(double[] data, double dataMin, double dataMax, int max_itr) {//k���ϖ@�B�e�f�[�^�_���ŏ��ɑ�����N���X�^�̃C���f�b�N�X��l�Ƃ��Ď��z���Ԃ��B
		int size = data.length;
		int[] cluster = new int[size];
		
		//�N���X�^���S
		double[] centers = new double[clusterMax];
		for (int k=0; k<clusterMax; k++) {
			centers[k] = dataMin + (dataMax - dataMin)*Math.random();
		}
		
		//�C�e���[�^
		for (int itr=0; itr<max_itr; itr++) {
			//�ł����S�̋߂��N���X�^�Ƀf�[�^��z��
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
			
			//���S�Čv�Z
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