<!doctype html>
<html lang="en">
  <head>
    <!-- Required meta tags -->
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">

    <!-- Bootstrap CSS -->
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.1/css/bootstrap.min.css" integrity="sha384-WskhaSGFgHYWDcbwN70/dfYBj47jz9qbsMId/iRN3ewGhXQFZCSftd1LZCfmhktB" crossorigin="anonymous">
    <title>Report for "$library" from $current</title>
  </head>
  <body>
    <div class="container mt-4">
    	<div class="row">
    		<div class="col">
	    		<h1>Report for "$library" from $current</h1>
    		</div>
    	</div>
    </div>
    <div class="container">
    	<div class="row">
    		<div class="col">
	    		<h2>How many reads in each replica?!</h2>
    		</div>
    	</div>
    	<div class="row">
    		<div class="col">
    			$output_rpl_count
    		</div>
    	</div>
    </div>
    <div class="container">
    	<div class="row">
	    	<div class="col">
	    		<h2>Data with unique barcodes</h2>
	    	</div>
	    </div>
	    <div class="row">
	    	<div class="col">
	    		$data
	    	</div>
    	</div>
    	<div class="row">
    		<div class="col">
		    	<h2>Aligning mapping, normalization and <i>"count"</i> expression</h2>
    		</div>
    	</div>
    	<div class="row">
	    	<div class="col">
	    		<p>Align <strong>mapping</strong> data with each other and return common combinations in set (bc, seq)</p>
	    		<p>Compare <strong>normalization</strong> replicates with each other and return common barcodes, but with its own value in each replica</p>
	    		$map_norm_1_2_data
	    	</div>
	    	<div class="col">
	    		<p>For <strong>expression</strong> - don't compare replicates, instead, we count the reads for each barcode considering mutated barcodes</p>
	    		$expr_1_2_data
	    	</div>
    	</div>
    </div>
    <div class="container my-3">
    	<div class="row">
    		<div class="col">
	    		<h2>Mapped data</h2>
    		</div>
    	</div>
    	<div class="row">
    		<div class="col">
    			<p>Count of Barcodes are common for mapping &amp; normalization data</p>
    			$mapped_ratio_data
    		</div>
    		<div class="col">
    			<p>Count of Barcodes are common for mapping &amp; normalization  <i>VS</i> expression - for each expression replica</p>
    			$mapped_norm_expression_data
    		</div>
    	</div>
    </div>
    <div class="container">
	    <div class="row">
	    	<div class="col">
		    	<h2>Control wt &amp; &Delta;C data</h2>
	    	</div>
	    </div>
	    <div class="row">
		    <div class="col">
			    $output_control_data
		    </div>
	    </div>
    </div>
    <div class="mb-4"></div>
    <div class="container mx-auto mb-3" style="width: 300px;">
	    <blockquote class="blockquote text-center">
			<footer class="blockquote-footer">Template developed by <cite title="Source Title">A.V. Ivankin</cite></footer>
		</blockquote>
	</div>
    <!-- Optional JavaScript -->
    <!-- jQuery first, then Popper.js, then Bootstrap JS -->
    <script src="https://code.jquery.com/jquery-3.3.1.slim.min.js" integrity="sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo" crossorigin="anonymous"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.3/umd/popper.min.js" integrity="sha384-ZMP7rVo3mIykV+2+9J3UJ46jBk0WLaUAdn689aCwoqbBJiSnjAK/l8WvCWPIPm49" crossorigin="anonymous"></script>
    <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.1.1/js/bootstrap.min.js" integrity="sha384-smHYKdLADwkXOn1EmN1qk/HfnUcbVRZyYmZ4qpPea6sjB/pTJ0euyQp0Mk8ck+5T" crossorigin="anonymous"></script>
  </body>
</html>