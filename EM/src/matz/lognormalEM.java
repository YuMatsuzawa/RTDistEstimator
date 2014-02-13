package matz;

/* EM�A���S���Y���Ń��f�����肷��N���X�B
 * �����LognormalMM�i�����ΐ����K���z�j�BclusterMax�͖{���I�ɂ̓p�����[�^�Ȃ̂ł������񂵂āA
 * ������܂߂ď��ʊ�������ōœK���f�������߂Ȃ��Ƃ����Ȃ��͂��B
 * �v�Z�����Ƃ���X�V�����\��ʂ�GMM�ƉZ��������̂ŃR�[�h�����s���Ԃ��Z��B
 * �Ƃ肠����RTCount�̏o�͕��@�ɍ��킹���̂ŁAretweetXXXXXXXXXXXX�f�B���N�g������XXXXXXXXXXXX.csv�t�@�C�������ׂēǂݍ��ނ͂��B
 * ����LNMM�ɂ͎��Ԓx��p�����[�^�͂܂���������Ă��Ȃ��B
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
	
	public static FilenameFilter regex(String regex){ //���K�\���ɂ��t�@�C�����t�B���^�[��Ԃ�
		final String regex_ = regex;
		return new FilenameFilter() {
			public boolean accept(File file, String name){
				boolean ret = name.matches(regex_);
				return ret;
			}
		};
	}
	
	//���s���Ԍv��
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
			
			//����Ƀf�[�^��T��
			File rootPath = new File(infroot);
			File[] dataDirs = rootPath.listFiles(regex("retweetLog[\\d]+"));
			for (File dir : dataDirs) {
				String[] dataFileNames = dir.list(regex("[\\d]+\\.csv"));
				for (String fileName : dataFileNames) {
					int resFileNum = 1;
					String resFileName = fileName.replaceAll("([\\d]+\\.)(csv)", "$1res\\.ln\\." + resFileNum + "\\.$2");
					File resDir = new File(dir, fileName.replaceAll("([\\d]+\\.)(csv)", "$1res"));
					if(!resDir.isDirectory()) resDir.mkdir();
					File resFile = new File(resDir, resFileName); //���ʏo�͐�
					if (clusteringMethod == KMEANS){ //�ő�N���X�^�����w�肵�āA�N���X�^���𑝂₵�Ȃ��玎�s����iK���ϖ@�ŃN���X�^�����j
						for (int cluster=1; cluster <= clusterMax; cluster++){
							while (resFile.exists()){
								//�o�̓t�@�C�������݂���(����1��ȏ㏈���ς�)�Ȃ�ԍ����C���N�������g����B
								resFileNum++;
								resFileName = fileName.replaceAll("([\\d]+\\.)(csv)", "$1res\\.ln\\." + resFileNum + "\\.$2");
								resFile = new File(resDir, resFileName);
							}
							if (!resFile.exists()) { 
								File dataFile = new File(dir, fileName); //���͌�
								new lognormalEM(cluster, dataFile, resFile);
							}
						}
					} else if(clusteringMethod == KLEINBERG) { //�o�[�X�g���m���g���ăN���X�^���𔻒肷��
						while (resFile.exists()){
							//�o�̓t�@�C�������݂���(����1��ȏ㏈���ς�)�Ȃ�ԍ����C���N�������g����B
							resFileNum++;
							resFileName = fileName.replaceAll("([\\d]+\\.)(csv)", "$1res\\.ln\\." + resFileNum + "\\.$2");
							resFile = new File(resDir, resFileName);
						}
						if (!resFile.exists()){
							File dataFile = new File(dir, fileName); //���͌�
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
			//�o�̓t�@�C�������݂���(����1��ȏ㏈���ς�)�Ȃ�ԍ����C���N�������g����B
			iclogFileNum++;
			iclogFileRenameTo = new File(resFile.getParentFile(), resFile.getName().replaceAll("([\\d]+\\.)(res\\.ln\\.[\\d]+\\.)(csv)", "$1" + icString +"\\." + iclogFileNum + "\\.$3"));
		}
		if(iclogFile.length() > 0) {
			iclogFile.renameTo(iclogFileRenameTo);
			stopWatch(iclogFileRenameTo.toString());
		}
		
	}
	
	public lognormalEM(String[] args) { //�������R���X�g���N�^
		try {
			/* �����̏����͌Œ�Ƃ���D
			 * �f�t�H���g����ύX���Ȃ��㑱�̈����͏ȗ����Ă��������C
			 * �O���̈������ȗ����Č���̈�����������͂��邱�Ƃ͂ł��Ȃ��D
			 * ����̈���������ύX�������ꍇ�͑O���̈����̃f�t�H���g�l��\�L����D
			 * �f�t�H���g��"./ 2 10000 0 1"
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
	
	public lognormalEM(int cluster, File dataFile, File resFile) throws FileNotFoundException, IOException { //�v�Z�R���X�g���N�^�B���Ă̒ʂ���̓f�[�^�t�@�C���Əo�͐�t�@�C�����w�肷��B
		List<Double> dataList = new ArrayList<Double>();
		
		stopWatch("Reading: " + dataFile.toString());
		/* RT���O�̕���1�s�ڂɃI���W�i��Tweet��URL��������Ă���̂ŁA�f�[�^��2�s�ڂ������B
		 * ����Ă��̃I�t�Z�b�g���w�肵�Ă������ƁB
		 */

		InputStreamReader isr = new InputStreamReader(new FileInputStream(dataFile));
		BufferedReader br = new BufferedReader(isr);
		String head = br.readLine(); //�擪�s
		
		String line = new String();
		while((line = br.readLine()) != null) {//�f�[�^�����X�g�ɓǂݍ��݁Bnull�i�t�@�C�������j�Ńu���C�N����B
			try {
				StringTokenizer st = new StringTokenizer(line, ",");
				dataList.add(Double.parseDouble(st.nextToken()));
			} catch(NumberFormatException e) { //Double�Ƀp�[�X�ł��Ȃ��ςȒl�����������烊�^�[�����Ă����B
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
		if (size < 100) { //1�T�Ԉȓ���RT����100�����Ȃ甲����B�኱���ꂷ����
			stopWatch("Disqualified. Aborting.");
			return;
		}
		//double[] data = new double[dataList.size()];
		double[] data = new double[size];
		int dataMin = Integer.MAX_VALUE;
		int dataMax = 0;
		int dataIndex = 0;
		for(int i=0; i<dataList.size(); i++){ //�����Ń��X�g��z��Ɉڂ��ւ��A�ȍ~�͂��̔z�񂩂�f�[�^���Q�Ƃ���Bdouble���v���~�e�B�u�^�Ȃ̂ł��傤�ǂ������\�b�h���Ȃ��炵��
			if (dataList.get(i) <= upperMax) {//�g�p����f�[�^���Ɍ��ōi��
				data[dataIndex] = dataList.get(i);
				//data[dataIndex] = dataList.get(i) + lowerOffset;
				if (data[dataIndex] < dataMin) dataMin = (int)data[dataIndex];
				if (data[dataIndex] > dataMax) dataMax = (int)data[dataIndex];
				dataIndex++;
			}
		}
		dataList.clear();
		
		//���������N���X�^������
		int[] initCluster = new int[size];
		if(cluster > 0) initCluster = kmeans(cluster, data, dataMin, dataMax, maxIter);
		else if(cluster == 0) {
			initCluster = kleinberg(cluster, data, dataMin, dataMax, maxIter);
			cluster = kleinbergCluster;
		}
		stopWatch("Processing EM for " + cluster + " cluster..");
		
		//�p�����[�^�̏�����
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
		
		//�C�e���[�^
		double[][] eAff = new double[size][cluster];
		final double cnst = Math.sqrt(2.0*Math.PI);
		ArrayList<Double> likelihoodSeries = new ArrayList<Double>();
		for(int itr=0; itr<maxIter; itr++){
			
			//E�X�e�b�v
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
					} catch (Exception e) { //�s�K�؂ȕ��U�����N���X�^��catch
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
			
			//M�X�e�b�v
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
			
			//�ޓx
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
		//���ʏo��
		OutputStreamWriter osw = new OutputStreamWriter(new FileOutputStream(resFile));
		BufferedWriter bw = new BufferedWriter(osw);
		bw.write(head);
		bw.newLine();
		
		//�ߎ��Ȑ��v���b�g�����E�o��
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
		
		//�N���X�^�����o��
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
		
		//�ޓx���ڏo��
		for(int i=0; i<likelihoodSeries.size(); i++){
			bw.write(likelihoodSeries.get(i).toString());
			bw.newLine();
		}
		bw.close();
		osw.close();
		
		//criterion�o��
		if (criterion == AIC) AIC(cluster, size, data, mean, var, weight, resFile);
		else if (criterion == BIC) BIC(cluster, size, data, mean, var, weight, resFile);
	}
	
	private void BIC(int cluster, int size, double[] data, double[] mean,
			double[] var, double[] weight, File resFile) {//BIC�o��
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
			int numOfParam = 3; //�N���X�^������̃p�����[�^��
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
			double[] var, double[] weight, File resFile) { //AIC�o��
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
			int numOfParam = 3; //�N���X�^������̃p�����[�^��
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

	public int[] kmeans(int cluster, double[] data, double dataMin, double dataMax, int max_itr) {//k���ϖ@�B�e�f�[�^�_���ŏ��ɑ�����N���X�^�̃C���f�b�N�X��l�Ƃ��Ď��z���Ԃ��B
		stopWatch("Clustering by K-means method..");
		
		int size = data.length;
		int[] initCluster = new int[size];
		
		//�N���X�^���S
		double[] centers = new double[cluster];
		for (int k=0; k<cluster; k++) {
			centers[k] = dataMin + k*(dataMax - dataMin)/cluster; //�����W���ɓ��Ԋu�ɔz�u
			//centers[k] = dataMin + (dataMax - dataMin)*Math.random(); //�����W���Ń����_���ɔz�u
		}
		
		//�C�e���[�^
		for (int itr=0; itr<max_itr; itr++) {
			int cutoff = size/100;
			
			int[] dataNum2 = new int[cluster]; 
			/* �ŏ��Ƀ����_���������ꂽ�ہA�ʒu�������āi�W�c���牓���āj�����o����l�ȉ���������U���Ȃ������N���X�^�����̂܂܂ɂȂ��Ă���ƁA��Ƀo�O�������N�����i���̃N���X�^�������U��0�ɂȂ��Ă��܂��̂Łj�B
			 * ���������邽�߂ɂ��������o����l�ȉ��̃N���X�^�����������炻�̃N���X�^�������S�𜓈ӓI�ɒ�������B
			 */
			//�ł����S�̋߂��N���X�^�Ƀf�[�^��z��
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
			
			//���S�Čv�Z
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
				if(dataNum[k]<=cutoff) continue; //�������ڂ�͌�ŋ~�ς���B�������ڂ�臒lcutoff�𓱓�
				centers[k] /= dataNum[k];
			}
			for(int k=0; k<cluster; k++) {//�������ڂ�~��
				double delta = 0.001 * (dataMax - dataMin);
				if(dataNum[k]<=cutoff) {
					int l = k;
					while (l==k || dataNum[l] <= cutoff) l = (int)(Math.random()*cluster);
					centers[k] = centers[l] + 0.5*delta;
					centers[l] = centers[l] - 0.5*delta; //���Ƀ����o�������Ă���N���X�^�������_���Ɏ���Ă��āA���̃N���X�^�̒��S�𒆐S�ɕ�delta�������ꂽ�ʒu�ɒ��S����点�Ă��炤���Ƃŋ~�ς���
				}
			}
		}
		return initCluster;
	}
	
	public int[] kleinberg(int cluster, double[] data, double dataMin, double dataMax, int max_itr){ //Kleinberg�̃o�[�X�g���m�A���S���Y���ɂ��N���X�^�����O
		stopWatch("Clustering by Kleinberg Algorithm..");
		
		int size = data.length;
		
		//frequency map
		double range = dataMax - dataMin;
		//int pitch = 10; //�p�x���z�̋�ԁi���j
		int pitch = 10*60; //�p�x���z�̋�ԁi�b�j
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
			//��ԃR�X�g
			double[] candidateCost = new double[stateMax];
			for(int s=0; s<stateMax; s++){
				candidateCost[s] = - Math.log(Math.pow(stateCost[s], freqInt) * Math.pow(1 - stateCost[s], denomInt - freqInt) * combination(denomInt, freqInt))
						+ sequentialCost;
				if(a>1 && state[a-1] < s) candidateCost[s] += (s - state[a-1]) * gamma; //��ԑJ�ڃR�X�g
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



