from flask import Flask, render_template, request, jsonify
from flask_cors import CORS
import os
from stl import mesh
import numpy as np

app = Flask(__name__)
CORS(app)
app.config['UPLOAD_FOLDER'] = 'uploads'
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/upload', methods=['POST'])
def upload_file():
    if 'stlFile' not in request.files:
        return jsonify({'error': 'No file uploaded'}), 400

    file = request.files['stlFile']
    filename = file.filename
    filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    file.save(filepath)

    # Load STL mesh
    model = mesh.Mesh.from_file(filepath)

    # Calculate surface area
    surface_area = np.sum(model.areas)

    # Volume
    volume, cog, inertia = model.get_mass_properties()

    # Bounding box (dimensions)
    minx, maxx = np.min(model.x), np.max(model.x)
    miny, maxy = np.min(model.y), np.max(model.y)
    minz, maxz = np.min(model.z), np.max(model.z)

    dimensions = {
        "length": round(maxx - minx, 4),
        "width": round(maxy - miny, 4),
        "height": round(maxz - minz, 4)
    }

    result = {
        "filename": filename,
        "surface_area_m2": round(surface_area, 6),
        "volume_m3": round(volume, 6),
        "dimensions_m": dimensions
    }

    return jsonify(result)


if __name__ == '__main__':
    app.run(host='127.0.0.1', port=5000, debug=True)