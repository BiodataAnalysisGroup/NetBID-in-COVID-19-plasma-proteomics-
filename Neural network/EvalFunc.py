import torch

from sklearn.metrics import roc_auc_score, f1_score
import sklearn.metrics as metrics

import matplotlib.pyplot as plt

def auc(y_true, y_pred):
	###if gpu is being used, transferring back to cpu
	if torch.cuda.is_available():
		y_true = y_true.cpu().detach()
		y_pred = y_pred.cpu().detach()
	###
	auc = roc_auc_score(y_true.detach().numpy(), y_pred.detach().numpy())



	fpr, tpr, threshold = metrics.roc_curve(y_true.detach().numpy().ravel(), y_pred.detach().numpy().ravel())
	roc_auc = metrics.auc(fpr, tpr)
	plt.figure()
	plt.title('Receiver Operating Characteristic')
	plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % auc)
	plt.legend(loc = 'lower right')
	plt.plot([0, 1], [0, 1],'r--')
	plt.xlim([0, 1])
	plt.ylim([0, 1])
	plt.ylabel('True Positive Rate')
	plt.xlabel('False Positive Rate')
	plt.savefig('auc.svg', metadata=None)

	return(auc)

def f1(y_true, y_pred):
	###covert one-hot encoding into integer
	y = torch.argmax(y_true, dim = 1)
	###estimated targets (either 0 or 1)
	pred = torch.argmax(y_pred, dim = 1)
	###if gpu is being used, transferring back to cpu
	if torch.cuda.is_available():
		y = y.cpu().detach()
		pred = pred.cpu().detach()
	###
	f1 = f1_score(y.detach().numpy(), pred.detach().numpy())
	return(f1)
